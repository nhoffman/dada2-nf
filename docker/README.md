# A Docker image providing dada2 (among other things)

This repository contains a Dockerfile for building a Docker image
inheriting from r-base with a specified release version of the dada2
package installed.

## Building

Create the Docker image locally like this:

```
./build.sh
```

An optional first argument can be used to specify the dada2 version.

## Version numbering

Images are tagged using the value of ``git describe --tags --dirty``
for this repository. The repo will have annotated tags corresponding
to the dada2 release version. A tagged image will therefore have the format

```
<dada2-version>[-<commits-since-tag>-<short-sha>]
```

For example:

```
dada2-nf:v1.12-5-g385b439
```

## Publishing images

Images are hosted at https://quay.io/nhoffman/dada2-nf

To push the most recent image:

```
image=$(docker images dada2-nf --format "{{.Repository}}:{{.Tag}}" | head -n1)
docker tag "$image" "quay.io/nhoffman/$image"
docker login quay.io/nhoffman
docker push "quay.io/nhoffman/$image"
```

