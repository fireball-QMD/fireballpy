.. currentmodule:: fireballpy

.. _fdata:

What are FDatas?
================

In the previous section we ended by noting a download
that happened as we tried to compute the energies
for our example.
This download is requiered as we need what we call *FDatas*.

So, what is a *FData*? It is the collection of precomputed
integrals of the base functions which allow for massive
speed improvement as they do not need to be calculated every time.
In exchange, it does occupy space in the hard disk (few hundred MBs)
and requiere a download the first time they are used.
You can check what *FDatas* are available by using the
``available_fdatas()`` function.

Before moving further we should address the question:
"Are you downloading files into my computer?" **Yes**.
We aknowledge the posible security concerns which is why
we perform sha256sum checks to ensure nothing wrong
happens with the download.
We also ask before updating any *FData* and replacing files.
Nonetheless, this a good time as any to remind that no
script (ours or not) should NEVER be executed with administrator
priviledges unless you place the utmost trust on it and it is strictly necessary.
To end with this discourse, you may check where these files
are downloaded with the ``get_fb_home()`` function.
The download path may be controlled with the environmental variable
``FIREBALL_HOME``. By default it will take the path
given by the environmental variable ``XDG_CACHE_HOME`` or
fallback to the Unix default ``~/.cache/``.

The next question one may have is what *FData* to use.
We have prepared a whole set of *FDatas* available
to be fetched from FireballPy.
Find a brief summary of them in the next table:

.. TODO table with FDatas

If you find that no *FData* may suit your need please feel free to contact us!

After this gentle introduction to the concept of *FDatas* is time to return
to coding and explore the nice things FireballPy can do!
