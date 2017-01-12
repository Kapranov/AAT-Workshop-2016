# AAT Observational Techniques Workshop 2016
### Data Archives Notes.

What is a data archive? - Vizier

## Vizier

Interested in the paper Green et al., [MNRAS 437, 1070
(2014)](http://adsabs.harvard.edu/abs/2014MNRAS.437.1070G).

Would be interested in plotting Halpha vs Mgas to test the analysis in
the paper.

Type numbers into python from the paper? \*Ridiculous!!!\*

## Use Vizier to retrieve tabular data for published papers

* Go to Vizier!
* Search "Green 2014"

Catalog appears (name is J/MNRAS/437/1070) in [search results][2].

Looking at the catalogs available, it looks like we want the "sample"
catalog.

We could download this file in some format, and then parse it into
python...or we could just see if python can get the data directly.

## Connecting to Vizier using a ``script/program``

Search "Vizier API" on Google. Looks like there is a python [package
called "Astroquery"][4] that can directly query Vizier.

* Install [Anaconda 4.2.0][5]

```bash
$ bash Anaconda2-4.2.0-Linux-x86_64.sh
```

* Configure [Conda][6] and Install AstroConda

```bash
$ conda config --add channels http://ssb.stsci.edu/astroconda
# Writes changes to ~/.condarc

$ conda create -n astroconda stsci

$ source activate astroconda
$ pip search Astroquery
$ pip install Astroquery
$ source deactivate astroconda
```

[1]: https://archive.gemini.edu
[2]:
http://vizier.u-strasbg.fr/viz-bin/VizieR-2?-ref=VIZ5727fab821b6&-to=2&-from=-2&-this=-2&%2F%2Fsource=&-out.max=50&%2F%2FCDSportal=&-out.form=HTML+Table&-out.add=_r&-out.add=_RAJ%2C_DEJ&%2F%2Foutaddvalue=&-sort=_r&-order=I&-oc.form=sexa&-meta.foot=1&-meta=1&-meta.ucd=2&-source=Green+2014&%21-2%3B=+Find...+&-ucd=&%2F%2Fucdform=on&-c=&-c.eq=J2000&-c.r=++2&-c.u=arcmin&-c.geom=r&-sort=_r&-order=I&-sort=_r&-order=I&-meta.ucd=2&-usenav=1&-bmark=GET
[3]:
https://archive.gemini.edu/searchform/cols=CTOWEQ/notengineering/NIFS/ra=22:17:39.85/dec=+00:15:26.42/NotFail
[4]: http://astroquery.readthedocs.io/en/latest/vizier/vizier.html

[5]: https://www.continuum.io/downloads
[6]:
http://astroconda.readthedocs.io/en/latest/installation.html#obtain-anaconda
