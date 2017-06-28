FROM alpine

MAINTAINER Philipp Rescheneder <philipp.rescheneder@gmail.com>

ARG VERSION_ARG

# Get most recent NextGenMap version from github
RUN apk --update upgrade && \
    apk add build-base gcc abuild binutils binutils-doc gcc-doc zlib-dev git cmake curl ca-certificates && \
    update-ca-certificates && \
    git clone https://github.com/Cibiv/NextGenMap.git && cd NextGenMap && git checkout $VERSION_ARG && mkdir -p build && cd build && cmake .. && make && cp ../bin/ngm-*/ngm /bin/ && cd .. && rm -rf NextGenMap && \
    apk del build-base gcc abuild binutils binutils-doc gcc-doc zlib-dev git cmake curl ca-certificates && \
    rm -rf /var/cache/apk/*
