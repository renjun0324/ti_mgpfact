FROM dynverse/dynwrap_latest:latest

ARG GITHUB_PAT

ENV JULIA_VERSION 1.6.6
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/$(echo $JULIA_VERSION | cut -d. -f 1-2)/julia-$JULIA_VERSION-linux-x86_64.tar.gz \
    && tar -xvzf julia-$JULIA_VERSION-linux-x86_64.tar.gz -C /usr/local/ \
    && ln -s /usr/local/julia-$JULIA_VERSION/bin/julia /usr/local/bin/julia \
    && rm julia-$JULIA_VERSION-linux-x86_64.tar.gz

RUN julia -e 'import Pkg; Pkg.add(url="https://github.com/renjun0324/MGPfact.jl")'

RUN julia -e 'import Pkg; Pkg.add(["Mamba", "RData", "JLD2", "Distributions"])'

RUN julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(name="RCall", version="0.13.15"))'

Run R -e 'devtools::install_version("JuliaCall","0.16")'

RUN R -e 'devtools::install_github("renjun0324/MURP@v0.6.5")'

RUN R -e 'devtools::install_cran(c("dplyr", "purrr", "stringr","JuliaCall", "pbmcapply", "doParallel", "reshape", "reshape2", "igraph", "graphlayouts","oaqc","parallelDist"))'

ARG CACHEBUST=12

COPY definition.yml run.R function.R /code/

RUN chmod +x /code/run.R

ENTRYPOINT ["/code/run.R"]

