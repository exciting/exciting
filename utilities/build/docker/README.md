# Docker

Docker is an open platform for developing, shipping, and running applications.
Docker enables you to separate your applications from your infrastructure so
you can deliver software quickly. Lots more information is available on the
[docker website](https://www.docker.com).

exciting uses docker to produce consistent environments for building and
testing in its continuous integration pipeline. Dockerfiles used to
generate the containers that our CI builds/tests run in, are provided here.

To build an image from a Dockerfile, install docker then type the following in
the docker (this) directory:

```bash
cp Dockerfile_* Dockerfile
docker build -t repository:tag . 
```

where `Dockerfile_*` is replaced with the appropriate filename. Repository is
replaced with the local repository name or git URL. For local use, one could
choose *exciting*. Finally, `tag` is replaced with any appropriate tag, for
example `debianBuster_GCC8`.


## Documentation Image

Dockerfile_docs provides the base image that contains all the dependencies
required to build exciting's documentation on the debian platform.

```bash
# Build an image containing all dependencies required to generate exciting’s
# documentation
docker build -t exciting:docs 
# Tag the image with the exciting_repository:name_of_image_in_repo
# The url used here is specific to exciting developers 
docker image tag exciting:docs gitdocker.physik.hu-berlin.de/sol/exciting:debian-base-docs
# Push the image such that it’s available to use in the CI 
docker push gitdocker.physik.hu-berlin.de/sol/exciting:debian-base-docs
```

## Debian (Buster) GCC8 Image

Dockerfile_debian-GCC8 provides all the base packages required to build and run exciting
with GCC 8, on Debian Buster (used in the SOL group at the time of writing). 

Alex notes that in the original docker image for GCC on to exciting's Gitlab (~ January 2021)
all tests pass with serial and MPISMP (unless marked as failing). This may not be the case
for Dockerfile_debian-GCC8, which was written in April 2022, and not fully-tested on the SMP version.
