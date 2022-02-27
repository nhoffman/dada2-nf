# A Docker image providing dada2 (among other things)

This repository contains a Dockerfile for building a Docker image
inheriting from r-base with a specified release version of the dada2
package installed.

## Building

Create the Docker image locally like this:

```
./build.sh
```

An optional first argument can be used to specify the dada2 version, and a
second can specify the git tag of this repository:

```
./build.sh 1.18 1.15
```

## Version numbering

Images are tagged using the value of ``git describe --tags --dirty``
for this repository. The repo will have annotated tags corresponding
to the dada2 release version. A tagged image will therefore have the format

```
<dada2-nf-tag>[-<commits-since-tag>-<short-sha>]
```

For example:

```
dada2-nf:v1.12-5-g385b439
```

## Publishing images

TODO: set up GitHub actions.

Images are hosted on Docker hub https://hub.docker.com/repository/docker/nghoffman/dada2-nf

To push the most recent image:

```
image=$(docker images dada2-nf --format "{{.Repository}}:{{.Tag}}" | head -n1)
docker tag "$image" "nghoffman/$image"
docker login nghoffman
docker push "nghoffman/$image"
```

