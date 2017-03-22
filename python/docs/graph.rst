.. _viz:

*************
Visualization
*************

This section describes how to visualize power networks using PFNET. To have this capability, PFNET needs the `Graphviz <http://www.graphviz.org/>`_ library ``libgvc``.

.. _viz_overview:

Overview
--------

To visualize a power network, a :class:`Graph <pfnet.Graph>` objects needs to be created. To do this, one needs to specify the power :class:`Network <pfnet.Network>` that is to be associated with the graph::

  >>> import pfnet

  >>> pfnet.ParserMAT().parse('ieee14.mat')

  >>> g = pf.Graph(net)

Then, a layout must be created for the graph. This can be done using the :class:`Graph <pfnet.Graph>` class method :func:`set_layout <pfnet.Graph.set_layout>`. This method uses the `sfdp algorithm of Graphviz <http://www.graphviz.org/content/root>`_. 

The :class:`Graph <pfnet.Graph>` class provides routines for coloring nodes (network buses) according to different criteria. For example, buses can be colored according to reactive power mismatches::

  >>> g.set_layout()

  >>> g.color_nodes_by_mismatch(pfnet.BUS_MIS_REACTIVE)

  >>> g.view()

.. image:: ./_static/graph.*
   :scale: 70%
   :align: center


