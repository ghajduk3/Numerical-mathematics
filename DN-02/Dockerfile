FROM julia:1.0

RUN apt-get -qq update && apt-get -y install git bzip2
WORKDIR /root
COPY Project.toml /root/ 
RUN julia --project=/root -e 'import Pkg; Pkg.instantiate()'
RUN mkdir -p /root/test/coverage
COPY test/coverage/Project.toml /root/test/coverage
RUN julia --project=/root/test/coverage -e 'import Pkg; Pkg.instantiate(); Pkg.build();'
RUN mkdir /root/docs
COPY docs/*.toml /root/docs/
RUN julia --project=/root/docs -e 'import Pkg; Pkg.instantiate(); Pkg.build();'