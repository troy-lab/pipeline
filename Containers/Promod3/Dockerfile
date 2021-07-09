FROM timm0049/ost:latest

# ARGUMENTS
###########
ARG PROMOD_VERSION="3.2.0"
ARG SRC_FOLDER="/usr/local/src"


# ENVIRONMENT
##############################################################################
ENV PROMOD_VERSION="${PROMOD_VERSION}"
ENV PROMOD_ROOT="/usr/local"


# INSTALL ProMod3
#################
RUN cd ${SRC_FOLDER} && \
    # copy promod release
    wget -O promod-${PROMOD_VERSION}.tar.gz -nc https://git.scicore.unibas.ch/schwede/ProMod3/-/archive/3.2.0/ProMod3-3.2.0.tar.gz && \
    mkdir promod-${PROMOD_VERSION} && \
    tar xf promod-${PROMOD_VERSION}.tar.gz -C ${SRC_FOLDER}/promod-${PROMOD_VERSION} --strip-components=1 && \
    mkdir -p ${SRC_FOLDER}/promod-${PROMOD_VERSION}/build && \
    cd ${SRC_FOLDER}/promod-${PROMOD_VERSION}/build && \
    # Build and install ProMod3
    cmake .. -DOST_ROOT=/usr/local \
             -DOPTIMIZE=1 \
             -DENABLE_SSE=1 \
             -DDISABLE_DOCUMENTATION=1 && \
 
    make && make check && make install && \
    # cleanup
    cd ${SRC_FOLDER} && rm ${SRC_FOLDER}/promod-${PROMOD_VERSION}.tar.gz && \
    rm -rf ${SRC_FOLDER}/promod-${PROMOD_VERSION}

ENTRYPOINT ["pm"]
