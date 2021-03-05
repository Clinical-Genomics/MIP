FROM python:3.8-slim

LABEL base_image="python/3.8-slim"
LABEL about.license="MIT License (MIT)"
LABEL about.version="c339d29"
LABEL about.home="https://github.com/northwestwitch/pdfmerger"
LABEL about.documentation="https://github.com/northwestwitch/pdfmerger"
LABEL about.tags="pdf,merge,bookmark,PyPDF2,reportlab,click,orientation,watermark"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN useradd --create-home --shell /bin/bash worker
WORKDIR /home/pdfmerger/app

RUN git clone https://github.com/northwestwitch/pdfmerger.git

WORKDIR /home/pdfmerger/app/pdfmerger

## Install
RUN git reset --hard c339d29 && \
    pip install -r requirements.txt && \
    pip install -e .

# Run commands as non-root user
RUN chown worker:worker -R /home/worker
RUN mkdir /home/pdfmerger/data
USER worker
ENTRYPOINT ["pdfmerger"]
CMD ["$@"]
