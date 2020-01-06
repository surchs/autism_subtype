# Based on jupyter notebook container
FROM jupyter/minimal-notebook

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

USER root

# ffmpeg for matplotlib anim
RUN apt-get update && \
    apt-get install -y --no-install-recommends ffmpeg && \
    rm -rf /var/lib/apt/lists/*

USER $NB_UID

RUN conda config --add channels conda-forge
RUN conda install -y ipywidgets \
                     jupyter \
                     numpy \
                     scipy \
                     scikit-learn \
                     matplotlib \
                     pandas \
                     nibabel \
                     statsmodels \
                     tqdm \
                     pillow \
                     nilearn && \
     conda remove --quiet --yes --force qt pyqt && \
     conda clean -tipsy && \
     # Activate ipywidgets extension in the environment that runs the notebook server
     jupyter nbextension enable --py widgetsnbextension --sys-prefix && \
     # Also activate ipywidgets extension for JupyterLab
     # Check this URL for most recent compatibilities
     # https://github.com/jupyter-widgets/ipywidgets/tree/master/packages/jupyterlab-manager
     jupyter labextension install @jupyter-widgets/jupyterlab-manager@^0.38.1 && \
     jupyter labextension install jupyterlab_bokeh@0.6.3 && \
     npm cache clean --force && \
     rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
     rm -rf /home/$NB_USER/.cache/yarn && \
     rm -rf /home/$NB_USER/.node-gyp && \
     fix-permissions $CONDA_DIR && \
     fix-permissions /home/$NB_USER

# Import matplotlib the first time to build the font cache.
ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot" && \
 fix-permissions /home/$NB_USER

USER $NB_UID