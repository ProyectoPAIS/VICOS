#!/bin/bash

#setup docker group based on hosts mount gid
echo "Adding hosts GID to docker system group"
# this only works if the docker group does not already exist
#DOCKER_SOCKET=/var/run/docker.sock
#DOCKER_GROUP=docker
#USER=docker
#
#DOCKER_GID=$(stat -c '%g' ${DOCKER_SOCKET})

#addgroup is distribution specific

#addgroup --gid ${DOCKER_GID} ${DOCKER_GROUP}
#adduser ${USER} ${DOCKER_GROUP}


exec "$@"