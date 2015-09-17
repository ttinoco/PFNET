.. _viz:

*************
Visualization
*************

This section describes how to visualize power networks using PFNET. To have this capability, PFNET needs the `Graphviz <http://www.graphviz.org/>`_ library.

.. _viz_overview:

Overview
--------

To visualize a power network, a ``Graph`` objects needs to be created. To do this, one needs to specify the power ``Network`` that is to be associated with the graph::

  >>> addpath(strcat(getenv('PFNET'),'/matlab'));

  >>> pfnet.load_library

  >>> net = pfnet.Network();
  >>> net.load('ieee14.mat');

  >>> g = pfnet.Graph(net);

Then, a layout must be created for the graph. This can be done using the ``Graph`` class method ``set_layout``. This method uses the `sfdp algorithm of Graphviz <http://www.graphviz.org/content/root>`_.  The graph can then be saved to a file in one of the `supported formats <http://www.graphviz.org/doc/info/output.html>`_ of `Graphviz <http://www.graphviz.org/>`_::

  >>> g.set_layout()

  >>> g.write('png','graph.png');

.. image:: ./_static/graph.*
   :scale: 70%
   :align: center


