Documentation
=============

This documentation is generated using `Sphinx <http://sphinx-doc.org/>`_ and
`Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ The two softwares are
bridged by means of the `Breathe extension <https://breathe.readthedocs.org/>`_
The online version of this documentation is built and served by `Read The Docs
<https://readthedocs.org/>`_.  The webpage http://mrchem.readthedocs.io/ is
updated on each push to the public GitHub repository.


How and what to document
------------------------

Doxygen enables documenting the code in the source code files thus removing a
"barrier" for developers.  To avoid that the code degenerates into a Big Ball
of Mud, it is mandatory to document directly within the source code classes and
functions.  To document general programming principles, design choices,
maintenance etc. you can create a .rst file in the ``doc`` directory. Remember
to refer the new file inside the ``index.rst`` file (it won't be parsed
otherwise).  Sphing uses `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ and `Markdown
<https://daringfireball.net/projects/markdown/>`_. Support for Markdown is not
as extensive as for reStructuredText, see `these comments
<https://blog.readthedocs.com/adding-markdown-support/>`_ Follow the guidelines
in :cite:`Wilson2014` regarding what to document.

Documeting methods in derived classes
-------------------------------------

Virtual methods should only be documented in the base classes.
This avoids unnecessary verbosity and conforms to the principle: "Document
_what_, not _how_" :cite:`Wilson2014`
If you feel the _how_ needs to be explicitly documented, add some notes in the
appropriate ``.rst`` file.

How does this work?
-------------------

To have an offline version of the documentation just issue ``make doc`` in the
build directory.  The HTML will be stored in ``doc/html``. Open the
``doc/html/index.html`` file with your browser to see and browse the
documentation.

.. warning::

   Building the documentation requires Doxygen, Sphinx, Perl and the Python
   modules PyYAML, Breathe and Matplotlib.
