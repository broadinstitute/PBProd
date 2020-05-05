FROM jupyter/datascience-notebook

RUN pip install --upgrade pip \
    && pip install datapackage \
    && pip install folium \
    && pip install --upgrade seaborn \
    && pip install geopandas \
    && pip install tensorflow \
    && pip install nbresuse \
    && pip install umap-learn

RUN jupyter serverextension enable --py nbresuse --sys-prefix
RUN jupyter nbextension install --py nbresuse --sys-prefix
RUN jupyter nbextension enable --py nbresuse --sys-prefix

RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager

RUN pip install scanpy
RUN pip install leidenalg louvain
RUN pip install umap-learn
RUN pip install anndata 
RUN pip install joblib 

RUN pip install opencv-python

#RUN pip install bokeh==1.4.0
RUN pip install bokeh
RUN jupyter labextension install @bokeh/jupyter_bokeh

RUN pip install "holoviews[recommended]"
RUN pip install graphviz
RUN pip install pysam

RUN pip install plotly==4.9.0
RUN jupyter labextension install jupyterlab-plotly@4.9.0
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager plotlywidget@4.9.0

USER root
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository universe
RUN apt-get update
RUN apt-get install -y graphviz
RUN apt-get install -y libgraphviz-dev
RUN apt-get install -y graphviz-dev 

RUN apt-get install -y graphviz libgraphviz-dev graphviz-dev pkg-config

RUN echo "jovyan:jovyan" | chpasswd

USER jovyan
RUN pip install pygraphviz
RUN jupyter labextension install @jupyterlab/toc

# Install SSW:
RUN mkdir ssw
RUN wget -O ssw-zip.zip https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/master.zip
RUN unzip ssw-zip.zip
RUN cd Complete-Striped-Smith-Waterman-Library-master/src && make libssw.so && mv libssw.so ../../ssw/. && mv ssw_lib.py ../../ssw/.
RUN touch ssw/__init__.py
RUN mv Complete-Striped-Smith-Waterman-Library-master ssw/.
RUN rm ssw-zip.zip 

# Install Tesserae2:
COPY tesserae2-master.zip /home/jovyan/
RUN unzip tesserae2-master.zip && cd tesserae2-master	&& pip install . && cd .. && rm tesserae2-master.zip

# Install parasail:
RUN pip install parasail

# Install abPOA:
USER root
RUN apt-get install -y libz-dev graphviz
USER jovyan
RUN pip install pyabpoa

# Install PyCytoscape
# (https://blog.jupyter.org/interactive-graph-visualization-in-jupyter-with-ipycytoscape-a8828a54ab63)
RUN pip install ipycytoscape
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-cytoscape

CMD ["start-notebook.sh"]

