# Readme

It appears that the github rendering of equations in Markdown is way too brittle to be useful. This is solved by generating a Jupyter notebook and displaying that instead. To do so, install Jupitext:

```
pip3 install jupytext
```

and then generate the Jupyter notebook:

```
jupytext --to notebook paper/paper.md
```
