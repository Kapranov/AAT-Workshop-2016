# AAT Observational Techniques Workshop 2016
### Data Archives Notes.

What is a data archive - Vizier?

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
astroquery (0.3.4)  - Functions and classes to access online data resources

$ pip install Astroquery
$ source deactivate astroconda
```

Start ``ipython`` and search catalog

```bash
$ conda info --env
$ source activate astroconda
(astroconda)$ ipython
```

```python
In [1]: from astroquery.vizier import Vizier

In [2]: cat_list = Vizier.get_catalogs('J/MNRAS/437/1070/sample')
In [3]: cat = cat_list[0]

In [4]: cat.colnames
In [5]: cat['logLIHa']
In [6]: cat['Mgas']

In [7]: import matplotlib.pyplot as plt
In [8]: plt.scatter(cat['Mgas'], cat['logLIHa'])
In [9]: plt.show()
```

Save image to ``figure_1.png``

## VizieR Queries (astroquery.vizier)

Astroquery is a set of tools for querying astronomical web forms and
databases. There are two other packages with complimentary functionality
as Astroquery: [astropy.vo][7]  is in the Astropy core and [pyvo][8] is
an Astropy affiliated package. They are more oriented to general

> Table Discover
> If you want to search for a set of tables, e.g. based on author name
> or other keywords, the [find_catalogs()][9] tool can be used:

```python
>>> from astroquery.vizier import Vizier
>>> catalog_list = Vizier.find_catalogs('Kang W51')
>>> print({k:v.description for k,v in catalog_list.items()})
{u'J/ApJ/706/83':   u'Embedded YSO candidates in W51 (Kang+, 2009)',
 u'J/ApJS/191/232': u'CO survey of W51 molecular cloud (Bieging+, 2010)'}
```

From this result, you could either get any of these as a complete
catalog or query them for individual objects or regions.

> Get a whole catalog
> If you know the name of the catalog you wish to retrieve, e.g. from
> doing a [find_catalogs()][9] search as above, you can then grab the
> complete contents of those catalogs:

```python
>>> catalogs = Vizier.get_catalogs(catalog_list.keys())
>>> print(catalogs)
TableList with 3 tables:
        '0:J/ApJ/706/83/ysos' with 33 column(s) and 50 row(s)
        '1:J/ApJS/191/232/table1' with 15 column(s) and 50 row(s)
        '2:J/ApJS/191/232/map' with 2 column(s) and 2 row(s)
```

Note that the row limit is set to 50 by default, so if you want to get a
truly complete catalog, you need to change that:

```python
>>> Vizier.ROW_LIMIT = -1
>>> catalogs = Vizier.get_catalogs(catalog_list.keys())
>>> print(catalogs)
>>> Vizier.ROW_LIMIT = 50
```

> Query an object
> For instance to query Sirius across all catalogs:

```python
>>> from astroquery.vizier import Vizier
>>> result = Vizier.query_object("sirius")
>>> print(result)
TableList with 275 tables:
        '0:METAobj' with 5 column(s) and 5 row(s)
        '1:ReadMeObj' with 5 column(s) and 5 row(s)
        '2:I/34/greenw2a' with 16 column(s) and 1 row(s)
        ...
```

All the results are returned as a [TableList][10] object. This is a
container for [Table][11] objects. It is basically an extension
to [OrderedDict][12] for storing a [Table][11] against its name.
To access an individual table from the [TableList][10] object:

```python
>>> interesting_table = result['IX/10A/cor_ros']
>>> print(interesting_table)

_RAJ2000 _DEJ2000   _r       _1RXS       Rank ... e_Count  roff       Com         IncDate   UpdDate
  deg      deg     arcm                       ...          arcm                   "Y:M:D"   "Y:M:D"
-------- -------- ----- ---------------- ---- ... -------- ---- ---------------- ---------- -------
101.2885 -16.7119 0.264 J064509.3-164241    2 ...     0.03  1.6           160058 1999-02-15
101.2854 -16.7174 0.127 J064509.3-164241   14 ...     0.02  1.3 UK200422FP_N1_P1 1999-02-15
101.3156 -16.7339 1.952 J064509.3-164241   20 ...    0.002  2.4           160058 1999-02-15

```

To do some common processing to all the tables in the returned
[TableList][10] object, do just what you would do for a python
dictionary:

```python
>>> for table_name in result.keys():
...     table = result[table_name]
...     # table is now an `astropy.table.Table` object
...     # some code to apply on table
```

> Query a region
> To query a region either the coordinates or the object name around
> which to query should be specified along with the value for the radius
> (or height/width for a box) of the region. For instance to query a
> large region around the quasar 3C 273:

```python
>>> from astroquery.vizier import Vizier
>>> from astropy.coordinates import Angle
>>> result = Vizier.query_region("3C 273", radius=Angle(0.1, "deg"), catalog='GSC')
>>> print(result)

TableList with 3 tables:
        '0:I/254/out' with 15 column(s) and 17 row(s)
        '1:I/271/out' with 23 column(s) and 53 row(s)
        '2:I/305/out' with 35 column(s) and 206 row(s)
```

Note that the radius may also be specified as a string in the format
expected by [Angle][13]. So the above query may also be written as:

```python
>>> result = Vizier.query_region("3C 273", radius="0d6m0s", catalog='GSC')
>>> print(result)

TableList with 3 tables:
        '0:I/254/out' with 15 column(s) and 17 row(s)
        '1:I/271/out' with 23 column(s) and 53 row(s)
        '2:I/305/out' with 35 column(s) and 206 row(s)
```

Or using angular units and quantities from [astropy.units][14]:

```python
>>> from astroquery.vizier import Vizier
>>> import astropy.units as u
>>> result = Vizier.query_region("3C 273", radius=0.1*u.deg, catalog='GSC')
>>> print(result)

TableList with 3 tables:
        '0:I/254/out' with 15 column(s) and 17 row(s)
        '1:I/271/out' with 23 column(s) and 53 row(s)
        '2:I/305/out' with 35 column(s) and 206 row(s)
```

As mentioned earlier, the region may also be mentioned by specifying the
height and width of a box. If only one of the height or width is
mentioned, then the region is treated to be a square having sides equal
to the specified dimension.

```python
>>> from astroquery.vizier import Vizier
>>> import astropy.units as u
>>> import astropy.coordinates as coord
>>> result = Vizier.query_region(coord.SkyCoord(ra=299.590, dec=35.201,
... unit=(u.deg, u.deg),
... frame='icrs'),
... width="30m",
... catalog=["NOMAD", "UCAC"])
>>> print(result)

TableList with 3 tables:
  '0:I/297/out' with 19 column(s) and 50 row(s)
  '1:I/289/out' with 13 column(s) and 50 row(s)
  '2:I/322A/out' with 24 column(s) and 50 row(s)
```

One more thing to note in the above example is that the coordinates may
be specified by using the appropriate coordinate object from
[astropy.coordinates][15]. Especially for ICRS coordinates, some support
also exists for directly passing a properly formatted string as the
coordinate. Finally the ``catalog`` keyword argument may be passed in
either [query_object()][16] or [query_region(][17] methods. This may be
a string (if only a single catalog) or a list of strings otherwise.

> Specifying keywords, output columns and constraints on columns
> To specify keywords on which to search as well as conditions on the
> output columns, an instance of the [VizierClass][18] class specifying
> these must be first created. All further queries may then be performed
> on this instance rather than on the Vizier class.

```python
v = Vizier(columns=['_RAJ2000', '_DEJ2000','B-V', 'Vmag', 'Plx'],
... column_filters={"Vmag":">10"}, keywords=["optical", "xry"])

WARNING: xry : No such keyword [astroquery.vizier.core]
```

Note that whenever an unknown keyword is specified, a warning is emitted
and that keyword is discarded from further consideration. The behavior
for searching with these keywords is the same as defined for the web
interface ([for details see here][19]). Now we call the different query
methods on this Vizier instance:

```python
>>> result = v.query_object("HD 226868", catalog=["NOMAD", "UCAC"])
>>> print(result)

TableList with 3 tables:
        '0:I/297/out' with 3 column(s) and 50 row(s)
        '1:I/289/out' with 3 column(s) and 18 row(s)
        '2:I/322A/out' with 3 column(s) and 10 row(s)

print(result['I/322A/out'])

_RAJ2000    _DEJ2000   Vmag
  deg         deg      mag
---------- ---------- ------
299.572419  35.194234 15.986
299.580291  35.176889 13.274
299.582571  35.185225 14.863
299.594172  35.179995 14.690
299.601402  35.198108 14.644
299.617669  35.186999 14.394
299.561498  35.201693 15.687
299.570217  35.225663 14.878
299.601081  35.233338 13.170
299.617995  35.205864 13.946
```

When specifying the columns of the query, sorting of the returned table
can be requested by adding ``+`` (or ``-`` for reverse sorting order) in
front of the column name. In the following example, the standard
(``"*"``) columns and the calculated distance column (``"_r"``) of the
2MASS catalog (II/246) are queried, 20 arcsec around HD 226868. The
result is sorted in increasing distance, as requested with the (``"+"``)
in front of ``"_r"``.

```python
>>> v = Vizier(columns=["*", "+_r"], catalog="II/246")
>>> result = v.query_region("HD 226868", radius="20s")
>>> print(result[0])

_r    _RAJ2000   _DEJ2000   RAJ2000    DEJ2000 ...
arcs     deg        deg        deg        deg  ...
------ ---------- ---------- ---------- ----------
 0.134 299.590280  35.201599 299.590280  35.201599
10.141 299.587491  35.203217 299.587491  35.203217
11.163 299.588599  35.198849 299.588599  35.198849
12.289 299.586356  35.200542 299.586356  35.200542
17.688 299.586254  35.197994 299.586254  35.197994
```

Note: The special column ``"*"`` requests just the default columns of a
catalog; ``"**"`` would request all the columns.

> Query with table
> A [Table][11] can also be used to specify the coordinates in a region
> query if it contains the columns ``_RAJ2000`` and ``_DEJ2000``. The
> following example starts by looking for AGNs in the Veron & Cety
> catalog with a ``Vmag`` between 10.0 and 11.0. Based on the result of
> this first query, guide stars with a ``Kmag`` brighter than 9.0 are
> looked for, with a separation between 2 and 30 arcsec. The column
> ``_q`` in the ``guide`` table is a 1-based index to the ``agn`` table
> (not the 0-based python convention).

```python
agn = Vizier(catalog="VII/258/vv10",
... columns=['*', '_RAJ2000', '_DEJ2000']).query_constraints(Vmag="10.0..11.0")[0]

print(agn)

_RAJ2000 _DEJ2000  Cl  nR      Name     ...
  deg      deg                          ...
-------- -------- --- --- -------------
 10.6846  41.2694   Q              M 31 ...
 60.2779 -16.1108   Q     NPM1G-16.0168 ...
 27.2387   5.9067   A   *      NGC  676 ...
 40.6696  -0.0131   A          NGC 1068 ...
139.7596  26.2697   A          NGC 2824 ...
147.5921  72.2792   A          NGC 2985 ...
173.1442  53.0678   A          NGC 3718 ...
184.9608  29.6139   A         UGC  7377 ...
185.0287  29.2808   A          NGC 4278 ...
186.4537  33.5467   A          NGC 4395 ...
192.7196  41.1194   A          NGC 4736 ...
208.3612  40.2831   A          NGC 5353 ...

>>> guide = Vizier(catalog="II/246",
      column_filters={"Kmag":"<9.0"}).query_region(agn,
        radius="30s", inner_radius="2s")[0]
>>> guide.pprint()


 _q  _RAJ2000   _DEJ2000    _r    RAJ2000
       deg        deg      arcs     deg
--- ---------- ---------- ------ ----------
  1  10.686015  41.269630  3.917  10.686015
  1  10.685657  41.269550  2.911  10.685657
  1  10.685837  41.270599  5.462  10.685837
  1  10.683263  41.267456  7.878  10.683263
  1  10.683465  41.269676  3.228  10.683465
  3  27.238636   5.906066  2.294  27.238636
  4  40.669277  -0.014225  4.214  40.669277
  4  40.668802  -0.013064  2.876  40.668802
  4  40.669219  -0.012236  3.399  40.669219
  4  40.670761  -0.012208  5.271  40.670761
  4  40.670177  -0.012830  2.293  40.670177
 11 192.721982  41.121040  8.751 192.721982
 11 192.721179  41.120201  5.163 192.721179
```

# Reference/API

## astroquery.vizier Package

### VizieR Query Tool

**Author**: Julien Woillez ([jwoillez@gmail.com](mailto:jwoillez@gmail.com))
This package is for querying the VizieR service, primarily hosted at:
[http://vizier.u-strasbg.fr][20]

Note: If the access to catalogues with VizieR was helpful for your
research work, the following acknowledgment would be appreciated:

```
This research has made use of the VizieR catalogue access tool, CDS,
Strasbourg, France.  The original description of the VizieR service was
published in A&AS 143, 23
```

### Classes

[Vizier Class][21]([columns, column_filters, ...])
[Conf][22] Configuration parameters for [astroquery.vizier][23].

### 10 Jan 2017 [Oleg G.Kapranov](mailto:lugatex@yahoo.com)

[1]: https://archive.gemini.edu
[2]: http://vizier.u-strasbg.fr/viz-bin/VizieR-2?-ref=VIZ5727fab821b6&-to=2&-from=-2&-this=-2&%2F%2Fsource=&-out.max=50&%2F%2FCDSportal=&-out.form=HTML+Table&-out.add=_r&-out.add=_RAJ%2C_DEJ&%2F%2Foutaddvalue=&-sort=_r&-order=I&-oc.form=sexa&-meta.foot=1&-meta=1&-meta.ucd=2&-source=Green+2014&%21-2%3B=+Find...+&-ucd=&%2F%2Fucdform=on&-c=&-c.eq=J2000&-c.r=++2&-c.u=arcmin&-c.geom=r&-sort=_r&-order=I&-sort=_r&-order=I&-meta.ucd=2&-usenav=1&-bmark=GET
[3]: https://archive.gemini.edu/searchform/cols=CTOWEQ/notengineering/NIFS/ra=22:17:39.85/dec=+00:15:26.42/NotFail
[4]: http://astroquery.readthedocs.io/en/latest/vizier/vizier.html
[5]: https://www.continuum.io/downloads
[6]: http://astroconda.readthedocs.io/en/latest/installation.html#obtain-anaconda
[7]: http://docs.astropy.org/en/latest/vo/index.html
[8]: https://pyvo.readthedocs.io/en/latest/
[9]: http://astroquery.readthedocs.io/en/latest/api/astroquery.vizier.VizierClass.html#astroquery.vizier.VizierClass.find_catalogs
[10]: http://astroquoery.readthedocs.io/en/latest/api/astroquery.utils.TableList.html#astroquery.utils.TableList
[11]: http://docs.astropy.org/en/latest/api/astropy.table.Table.html#astropy.table.Table
[12]: http://docs.python.org/2/library/collections.html#collections.OrderedDict
[13]: http://docs.astropy.org/en/latest/api/astropy.coordinates.Angle.html#astropy.coordinates.Angle
[14]: http://docs.astropy.org/en/latest/units/index.html#module-astropy.units
[15]: http://docs.astropy.org/en/latest/coordinates/index.html#module-astropy.coordinates
[16]: http://astroquery.readthedocs.io/en/latest/api/astroquery.vizier.VizierClass.html#astroquery.vizier.VizierClass.query_object
[17]: http://astroquery.readthedocs.io/en/latest/api/astroquery.vizier.VizierClass.html#astroquery.vizier.VizierClass.query_region
[18]: http://astroquery.readthedocs.io/en/latest/api/astroquery.vizier.VizierClass.html#astroquery.vizier.VizierClass
[19]: http://vizier.u-strasbg.fr/vizier/vizHelp/1.htx
[20]: http://vizier.u-strasbg.fr/
[21]: http://astroquery.readthedocs.io/en/latest/api/astroquery.vizier.VizierClass.html#astroquery.vizier.VizierClass
[22]: http://astroquery.readthedocs.io/en/latest/api/astroquery.vizier.Conf.html#astroquery.vizier.Conf
[23]: http://astroquery.readthedocs.io/en/latest/vizier/vizier.html#module-astroquery.vizier

