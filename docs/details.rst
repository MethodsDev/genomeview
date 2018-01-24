.. _details:

GenomeView Details
==================

In the :ref:`previous section <tutorial>`, we saw how the :py:func:`genomeview.visualize_data()` quickly assembled a document to visualize read data. Here, we'll discuss the behind-the-scenes work of setting up a document, views and adding tracks. These steps are useful if you wish to customize any aspect of the visualization process.

.. contents:: :local:


Step 1: creating the document
-----------------------------

First we'll need a document. The only argument is to define the width of the view (think of this as its pixel width)::
    
    doc = genomeview.Document(900)

The document is where all the genome views will end up.


Step 2: creating the genome views
---------------------------------

We're starting to get into the action here -- a genome view defines a set of coordinates to visualize, and allows the addition of a number of tracks displaying different types of data for those coordinates.

To create a genome view, you'll first create a genome "source" (basically a link to the reference genome sequence), then derive a view with the coordinates you'd like to visualize::
    
    source = genomeview.FastaGenomeSource("/path/to/hg19.fasta")
    view = genomeview.GenomeView("locus 1", "chr1", 219158937, 219169063, "+", source)
    doc.add_view(view)

Note that all views and tracks take as their first argument a label which must be unique within the containing document or view.

You can add as many genome views as you'd like to a single document, allowing you to visualize multiple genomic loci in the same document.


Step 3: adding the tracks to the genome view
--------------------------------------------

The next step is to create tracks visualizing the actual data and add them to the genome view. Tracks are visualized in the order that they're added, so if you'd like to put the axis at the top, add it first, and if you'd like it at the bottom, add it last.

For example::

    bam_track_hg002 = genomeview.SingleEndBAMTrack("HG002", "/path/to/hg002.sorted.bam")
    view.add_track(bam_track_hg002)

    axis_track = genomeview.Axis("axis")
    view.add_track(axis_track)


Step 4: exporting the visualization
-----------------------------------

As mentioned in the previous section, the document can easily be visualized in-line in jupyter.

In addition, documents can be exported PDF or PNG files using the :py:func:`genomeview.save()` by changing the file suffix to ".pdf"/".png".

Note that this conversion requires the installation of `inkscape <https://inkscape.org/>`_, `libRsvg <https://wiki.gnome.org/action/show/Projects/LibRsvg>`_ or `webkitToPDF <https://github.com/nspies/webkitToPDF>`_ are installed.