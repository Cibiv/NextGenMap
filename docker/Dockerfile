FROM ubuntu

MAINTAINER Philipp Rescheneder <philipp.rescheneder@gmail.com>

ARG VERSION_ARG

# Get most recent NextGenMap version from github
RUN buildDeps='git wget gcc g++ libc6-dev make cmake zlib1g-dev ca-certificates' && \
    apt-get update && apt-get install -y $buildDeps --no-install-recommends && \
    update-ca-certificates && \
    git clone https://github.com/Cibiv/NextGenMap.git && cd NextGenMap && git checkout $VERSION_ARG && mkdir -p build && cd build && cmake .. && make && cp -r ../bin/ngm-*/* /usr/bin/ && cd .. && rm -rf NextGenMap && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get purge -y --auto-remove $buildDeps
