FROM quay.io/podman/stable:v2.1.1

LABEL org.opencontainers.image.source="https://github.com/eriksjolund/slurm-container-cluster" \
      org.opencontainers.image.title="slurm-container-cluster" \
      org.opencontainers.image.description="Slurm container cluster with norouter on Fedora"

ARG SLURM_TAG=slurm-20-11-0-1
ARG GOSU_VERSION=1.12

RUN dnf -y update \
    && dnf -y install \
       wget \
       bzip2 \
       perl \
       cpanminus \
       gcc \
       gcc-c++\
       git \
       gnupg \
       make \
       munge \
       munge-devel \
       python3-pip \
       mariadb-server \
       mariadb-devel \
       psmisc \
       bash-completion \
       vim-enhanced \
       procps-ng \
       net-tools \
       hostname \
    && dnf clean all \
    && rm -rf /var/cache/yum \
    && curl -fsSL https://github.com/norouter/norouter/releases/download/v0.6.1/norouter-$(uname -s)-$(uname -m).tgz | sudo tar xzvC /usr/local/bin \
#    && curl -o /usr/local/bin/norouter --fail -L https://github.com/norouter/norouter/releases/download/v0.6.1/norouter-$(uname -s)-$(uname -m) \
    && chmod 755 /usr/local/bin/norouter

RUN pip3 install Cython nose

RUN set -x \
    && wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-amd64" \
    && wget -O /usr/local/bin/gosu.asc "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-amd64.asc" \
    && export GNUPGHOME="$(mktemp -d)" \
    && gpg --keyserver ha.pool.sks-keyservers.net --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4 \
    && gpg --batch --verify /usr/local/bin/gosu.asc /usr/local/bin/gosu \
    && rm -rf "${GNUPGHOME}" /usr/local/bin/gosu.asc \
    && chmod +x /usr/local/bin/gosu \
    && gosu nobody true

RUN set -x \
    && git clone https://github.com/SchedMD/slurm.git \
    && pushd slurm \
    && git checkout tags/$SLURM_TAG \
    && ./configure --enable-debug --prefix=/usr --sysconfdir=/etc/slurm \
        --with-mysql_config=/usr/bin  --libdir=/usr/lib64 \
    && make install \
    && install -D -m644 etc/cgroup.conf.example /etc/slurm/cgroup.conf.example \
    && install -D -m644 etc/slurm.conf.example /etc/slurm/slurm.conf.example \
    && install -D -m644 etc/slurmdbd.conf.example /etc/slurm/slurmdbd.conf.example \
    && install -D -m644 contribs/slurm_completion_help/slurm_completion.sh /etc/profile.d/slurm_completion.sh \
    && popd \
    && rm -rf slurm \
    && groupadd -r slurm \
    && useradd -r -g slurm slurm \
    && mkdir /etc/sysconfig/slurm \
        /var/spool/slurmd \
        /var/run/slurmd \
        /var/run/slurmdbd \
        /var/lib/slurmd \
        /var/log/slurm \
        /data \
    && touch /var/lib/slurmd/node_state \
        /var/lib/slurmd/front_end_state \
        /var/lib/slurmd/job_state \
        /var/lib/slurmd/resv_state \
        /var/lib/slurmd/trigger_state \
        /var/lib/slurmd/assoc_mgr_state \
        /var/lib/slurmd/assoc_usage \
        /var/lib/slurmd/qos_usage \
        /var/lib/slurmd/fed_mgr_state \
    && chown -R slurm:slurm /var/*/slurm* \
    && /sbin/create-munge-key \
    chown munge:munge /var/run/munge && \
    sed -i 's/socket=\/var\/lib\/mysql\/mysql.sock/socket=\/var\/run\/mysqld\/mysqld.sock/g' /etc/my.cnf

COPY slurm.conf /etc/slurm/slurm.conf
COPY slurmdbd.conf /etc/slurm/slurmdbd.conf

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod 755 /usr/local/bin/docker-entrypoint.sh

#ENV LC_ALL en_US.UTF-8
#ENV LANG en_US.UTF-8
#ENV LANGUAGE en_US.UTF-8

RUN cpanm install Carton

# Add MIP into image
ENV MIP_INSTALL_DIR=/workspace/bin
COPY . "$MIP_INSTALL_DIR"
WORKDIR "$MIP_INSTALL_DIR"

# Make executable and add to binary to PATH
RUN chmod a+x "$MIP_INSTALL_DIR/mip"
ENV PATH "$PATH:$MIP_INSTALL_DIR"

# Remove any outside local from previous carton installs
RUN rm -rf local

RUN carton install

ENV PERL5LIB="$MIP_INSTALL_DIR"/local/lib/perl5

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

CMD ["slurmdbd"]
