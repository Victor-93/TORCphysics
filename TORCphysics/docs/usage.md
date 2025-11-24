Usage
=====

Installation
------------

To use TORCphysics, first install it using pip:

```console
(.venv) $ pip install TORCphysics
```

Creating recipes
----------------

To retrieve a list of random ingredients, you can use the
`TORCphysics.Circuit()` function:

::: TORCphysics.Circuit
    options:
      show_root_heading: true

<br>

The `kind` parameter should be either `"meat"`, `"fish"`, or `"veggies"`.
Otherwise, [`get_random_ingredients`][TORCphysics.get_random_ingredients] will raise an exception [`TORCphysics.InvalidKindError`](/api#lumache.InvalidKindError).

For example:

```python
>>> import TORCphysics
>>> TORCphysics.Circuit()
['shells', 'gorgonzola', 'parsley']
```