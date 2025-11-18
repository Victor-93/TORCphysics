Usage
=====

.. _installation:

Installation
------------

To use TORCphysics, first install it using pip:

.. code-block:: console

   (.venv) $ pip install TORCphysics

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``TORCphysics.get_random_ingredients()`` function:

.. autofunction:: TORCphysics.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`TORCphysics.get_random_ingredients`
will raise an exception.

.. autoexception:: TORCphysics.InvalidKindError

For example:

>>> import TORCphysics
>>> TORCphysics.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

