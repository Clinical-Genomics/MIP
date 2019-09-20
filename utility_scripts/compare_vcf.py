import cyvcf2
import copy
import logging
import sys
from tabulate import tabulate
from itertools import chain, tee

"""
Script compare two vcf-files on a number of selected measures

Usage: compare_vcf.py <vcf_1> <vcf_2>

Prints aggregated statistics on the compared fields and detailed comparison
of each found common variant in the vcf-files
"""

__version__ = '1.0.0'

LOG = logging.getLogger(__name__)

class Genotype:
    """Make string representation of cyvcf2 genotype"""
    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]
    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)
    __repr__ = __str__


# Conversion functions
get_score = lambda score: int(score.split(':')[-1])
get_alt_DP = lambda DP: [int(dp) if int(dp) >= 0 else None for dp in list(DP)]
get_alt_AD = lambda DP: [[int(i) if int(i) >= 0 else None for i in dp][-1] for dp in list(DP)]
get_GQ = lambda GQ: [int(gq[0]) if int(gq[0]) >= 0 else None for gq in GQ]
get_genotypes = lambda genotypes: [str(Genotype(gt)) for gt in genotypes]
get_qual = lambda qual: float(qual) if isinstance(qual, float) else str(qual)


# INFO_keys dictionary states what fields in the vcf that should be compared
# To add a field/fields that should be compared between two VCFs, add
#   NAME: {'type_conv': function, 'ID': (ID1,ID2), 'column': 'INFO'|FORMAT}
# Where NAME is a chosen name for the value, function a function that converts
# the value from cyvcf2 to something else, IDs are the name of the fields in the
# VCFs, these may be more than one, if the IDs in the two VCFs are different for
# The same measure. and column name, which is either of 'INFO' or 'FORMAT'
COMPARE_FIELDS = {
    'CHROM': {'type_conv': str, 'ID': ('CHROM',)},
    'POS': {'type_conv': int, 'ID': ('POS',)},
    'END': {'type_conv': int, 'ID': ('END',), 'column': 'INFO'},
    'RANK_SCORE': {'type_conv': get_score, 'ID': ('MutaccRankScore', 'RankScore'), 'column': 'INFO'},
    'TYPE': {'type_conv': str, 'ID': ('TYPE', 'SVTYPE', 'type',), 'column': 'INFO'},
    'DEPTH': {'type_conv': get_alt_DP, 'ID': ('DP',), 'column': 'FORMAT'},
    'ALT DEPTH': {'type_conv': get_alt_AD, 'ID': ('AD', 'AF',), 'column': 'FORMAT'},
    'GQ': {'type_conv': get_GQ, 'ID': ('GQ',), 'column': 'FORMAT'},
    'GENOTYPES': {'type_conv': get_genotypes, 'ID': ('GT',)},
    'QUAL': {'type_conv': float, 'ID': ('QUAL',)},
    'SV_LEN': {'type_conv': int, 'ID': ('SVLEN',), 'column': 'INFO'}
}


class VCFHandler(cyvcf2.VCF):
    """ Class to handle cyvcf2.VCF objects """
    def intersection(self, vcf):
        """ Finds intersection, i.e. common variants with another VCFHandler object"""
        for variant in self:
            if variant.INFO.get('SVTYPE', None) is None:
                variant_id = make_record_id(variant)
                record = find_record(vcf, record_id=variant_id)
            else:
                record = find_similar_sv(vcf, variant)
            if record is not None:
                yield (variant, record)

    def compare_fields(self, vcf, vcf_keys):
        """Compares self with another VCFHandler object"""
        intersection = self.intersection(vcf)
        for variants in intersection:
            record_1 = variants[0]
            record_2 = variants[1]
            comp = compare_records(record_1, record_2, vcf_keys)
            variant_id = make_record_id(record_1)
            comp['id'] = variant_id
            yield comp


def find_record(vcf: VCFHandler, record_id=None):
    """Find record in VCFHandler object based on record_id"""
    record_fields = record_id.split('_')
    region_str = f'{str(record_fields[0])}:{str(record_fields[1])}-{str(int(record_fields[1])+1)}'
    for record in vcf(region_str):
        if record_id == make_record_id(record):
            return record
    return None


def find_similar_sv(vcf, variant_record):
    """Given a cyvcf2 record, find closest matching SV in specified VCFHandler object"""
    pos = variant_record.POS
    end = variant_record.INFO.get('END')
    if end is None:
        return None
    type_1 = variant_record.INFO.get('SVTYPE')
    if type_1 is None:
        return None
    search_interval = 3000
    closest_match = None
    for hit in vcf(f"{variant_record.CHROM}:{pos-search_interval}-{pos+search_interval}"):
        type_2 = hit.INFO.get('SVTYPE')
        if type_2 is None or hit.INFO.get('END') is None:
            continue
        if hit.POS < pos - search_interval or hit.POS > pos + search_interval:
            continue
        if type_1 in type_2 or type_2 in type_1:
            if closest_match is None:
                closest_match = hit
                continue
            else:
                this_distance = abs(hit.POS-pos) + abs(hit.INFO['END']-end)
                match_distance = abs(closest_match.POS-pos) + abs(closest_match.INFO['END']-end)
                if this_distance <= match_distance:
                    closest_match = hit
    return closest_match


def make_record_id(record):
    """Make record id based on CHROM POS REF ALT"""
    record_id = '{}_{}_{}_{}'.format(
        record.CHROM,
        record.POS,
        record.REF,
        record.ALT[-1]
    )
    return record_id


def compare_records(record_1, record_2, vcf_keys=None):
    """ Compare two vcf records based on the keys specified in vcf_keys"""
    comparison = dict()
    for key, value in vcf_keys.items():
        record_1_key = None
        record_2_key = None
        for alt_key in value['ID']:
            if value.get('column', None) is None:
                if alt_key == 'POS':
                    record_1_key = record_1.POS
                    record_2_key = record_2.POS
                if alt_key == 'QUAL':
                    record_1_key = record_1.QUAL
                    record_2_key = record_2.QUAL
                if alt_key == 'GT':
                    record_1_key = record_1.genotypes
                    record_2_key = record_2.genotypes
                if alt_key == 'CHROM':
                    record_1_key = record_1.CHROM
                    record_2_key = record_2.CHROM
            elif value['column'] == 'INFO':
                if record_1_key is None:
                    record_1_key = record_1.INFO.get(alt_key)
                if record_2_key is None:
                    record_2_key = record_2.INFO.get(alt_key)
            elif value['column'] == 'FORMAT':
                if record_1_key is None:
                    try:
                        record_1_key = record_1.format(alt_key)
                    except KeyError:
                        record_1_key = None
                if record_2_key is None:
                    try:
                        record_2_key = record_2.format(alt_key)
                    except KeyError:
                        record_2_key = None
        if (record_1_key is not None) and (record_2_key is not None):
            comparison[key] = (value['type_conv'](record_1_key), value['type_conv'](record_2_key))
        else:
            log_msg = f"key {' or '.join(value['ID'])} not found in one of the records"
            LOG.debug(log_msg)
    return comparison


def aggregate_stats(comparisons: list):
    """Get aggregated statistics on a list of comparisons"""

    def compare_gt_list(gt_1, gt_2):

        for sample_1, sample_2 in zip(gt_1, gt_2):
            unphased_1 = str(sample_1).replace('|', '/')
            unphased_2 = str(sample_2).replace('|', '/')
            if unphased_1 == unphased_2:
                continue
            elif '.' in unphased_1 or '.' in unphased_2:
                continue
            else:
                return False
        return True

    agg_dict = {}
    for comp in comparisons:
        for key, value in comp.items():
            if key not in agg_dict.keys():
                if isinstance(value[0], (float, int)):
                    agg_dict[key] = {'compared_no': 1,
                                     'tot_sum_diff': value[1] - value[0],
                                     'higher': 0,
                                     'lower': 0,
                                     'same': 0}
                    if value[1] > value[0]:
                        agg_dict[key]['higher'] += 1
                    elif value[1] < value[0]:
                        agg_dict[key]['lower'] += 1
                    elif value[1] == value[0]:
                        agg_dict[key]['same'] += 1

                if isinstance(value[0], str):
                    agg_dict[key] = {'compared_no': 1,
                                     'same': 0,
                                     'different': 0}
                    if value[0] == value[1]:
                        agg_dict[key]['same'] += 1
                    else:
                        agg_dict[key]['different'] += 1

                if isinstance(value[0], list):
                    agg_dict[key] = {'compared_no': 1,
                                     'same': 0,
                                     'different': 0}
                    if compare_gt_list(value[0], value[1]):
                        agg_dict[key]['same'] += 1
                    else:
                        agg_dict[key]['different'] += 1
            else:
                if isinstance(value[0], (float, int)):
                    agg_dict[key]['compared_no'] += 1
                    agg_dict[key]['tot_sum_diff'] += value[1] - value[0]
                    if value[1] > value[0]:
                        agg_dict[key]['higher'] += 1
                    elif value[1] < value[0]:
                        agg_dict[key]['lower'] += 1
                    elif value[1] == value[0]:
                        agg_dict[key]['same'] += 1

                if isinstance(value[0], str):
                    agg_dict[key]['compared_no'] += 1
                    if value[0] == value[1]:
                        agg_dict[key]['same'] += 1
                    else:
                        agg_dict[key]['different'] += 1

                if isinstance(value[0], list):
                    agg_dict[key]['compared_no'] += 1
                    if compare_gt_list(value[0], value[1]):
                        agg_dict[key]['same'] += 1
                    else:
                        agg_dict[key]['different'] += 1

    for field, values in agg_dict.items():
        if 'tot_sum_diff' in values.keys():
            agg_dict[field]['mean_diff'] = agg_dict[field]['tot_sum_diff']/agg_dict[field]['compared_no']
    return agg_dict


def print_agg(agg_dict):
    """ Print aggregated statistics"""
    for field, values in agg_dict.items():
        print(f"{field}")
        for key, value in values.items():
            print(f"\t{key}: {value}")


def main():
    """ main function """
    try:
        positive_controls = sys.argv[1]
        vcf_file = sys.argv[2]
    except IndexError:
        LOG.info('Please provide truth set and vcf_file')
        raise
    positive_controls = VCFHandler(positive_controls)
    vcf_file = VCFHandler(vcf_file)
    compared_variants = positive_controls.compare_fields(vcf=vcf_file, vcf_keys=COMPARE_FIELDS)
    count = 0
    compared_variants_1, compared_variants_2 = tee(compared_variants)
    agg_dict = aggregate_stats(compared_variants_1)
    print_agg(agg_dict)
    for i in compared_variants_2:
        print(i['id'])
        header_names = list(set(i.keys())-{'id'})
        summary1 = ['expected'] + [str(i[key][0]) for key in header_names]
        summary2 = ['actual'] + [str(i[key][1]) for key in header_names]
        print(tabulate([summary1, summary2], headers=['VCF'] + header_names))
        print('\n')
        count += 1
        if count > 1000:
            break
    positive_controls.close()
    vcf_file.close()

# main script
if __name__ == '__main__':

    main()
