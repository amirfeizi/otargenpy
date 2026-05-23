otargenpy
=========

**Tidy Python interface to the Open Targets Platform GraphQL API.**

Query genes, diseases, drugs, variants, and genetic evidence directly from Python
and receive analysis-ready pandas DataFrames.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   api


Installation
------------

From PyPI::

   pip install otargenpy

From GitHub (latest)::

   pip install git+https://github.com/amirfeizi/otargenpy.git

Quick Start
-----------

.. code-block:: python

   import otargenpy as ot

   # Adverse events for imatinib
   ae = ot.adverse_events_query("CHEMBL941")

   # Protein interactions for TP53
   inter = ot.interactions_query("ENSG00000141510", source_database="intact")

   # Plot results
   ot.plot_adverse_events(ae)
   ot.plot_interactions(inter)

Links
-----

- **GitHub:** https://github.com/amirfeizi/otargenpy
- **R sister package:** https://github.com/amirfeizi/otargen
- **Open Targets Platform:** https://platform.opentargets.org

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
