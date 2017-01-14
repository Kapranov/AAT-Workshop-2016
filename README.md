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

>>> guide = Vizier(catalog="II/246",
      column_filters={"Kmag":"<9.0"}).query_region(agn,
        radius="30s", inner_radius="2s")[0]
>>> guide.pprint()
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

# A Gallery of Queries

A series of queries folks have performed for research or for kicks.

## Example #1

This illustrates querying Vizier with specific keyword, and the use of
[astropy.coordinates][15] to describe a query. Vizier’s keywords can indicate
wavelength & object type, although only object type is shown here.

```python
>>> from astroquery.vizier import Vizier
>>> v = Vizier(keywords=['stars:white_dwarf'])
>>> from astropy import coordinates
>>> from astropy import units as u
>>> c = coordinates.SkyCoord(0,0,unit=('deg','deg'),frame='icrs')
>>> result = v.query_region(c, radius=2*u.deg)
>>> print len(result)
43
>>> result[0].pprint()
```

## Example #2

This illustrates adding new output fields to SIMBAD queries. Run
[list_votable_fields][24] to get the full list of valid fields.

```python
>>> from astroquery.simbad import Simbad
>>> s = Simbad()
>>> s.add_votable_fields('bibcodelist(2003-2013)')
>>> r = s.query_object('m31')
>>> r.pprint()
```
## Example #3

This illustrates finding the spectral type of some particular star.

```python
>>> from astroquery.simbad import Simbad
>>> customSimbad = Simbad()
>>> customSimbad.add_votable_fields('sptype')
>>> result = customSimbad.query_object('g her')
>>> result['MAIN_ID'][0]
'* g Her'
>>> result['SP_TYPE'][0]
'M6III'
```

## Example #4

```python
>>> from astroquery.simbad import Simbad
>>> customSimbad = Simbad()
>>> customSimbad.add_votable_fields('ra(d)','dec(d)')
>>> customSimbad.remove_votable_fields('coordinates')
>>> from astropy import coordinates
>>> C = coordinates.SkyCoord(0,0,unit=('deg','deg'), frame='icrs')
>>> result = customSimbad.query_region(C, radius='2 degrees')
>>> result[:5].pprint()
```

## Example #5

This illustrates a simple usage of the ``open_exoplanet_catalogue`` module.
Finding the mass of a specific planet:

```python
>>> from astroquery import open_exoplanet_catalogue as oec
>>> from astroquery.open_exoplanet_catalogue import findvalue
>>> cata = oec.get_catalogue()
>>> kepler68b = cata.find(".//planet[name='Kepler-68 b']")
>>> print findvalue( kepler68b, 'mass')
0.026 +0.008 -0.007
```

## Example #6

Grab some data from ALMA, then analyze it using the Spectral Cube
package after identifying some spectral lines in the data.

```python
from astroquery.alma import Alma
from astroquery.splatalogue import Splatalogue
from astroquery.simbad import Simbad
from astropy import units as u
from astropy import constants
from spectral_cube import SpectralCube

m83table = Alma.query_object('M83', public=True)
# Need knows of the fields - Member ous id
m83urls = Alma.stage_data(m83table['Member ous id'])
# Sometimes there can be duplicates: avoid them with
# list(set())
m83files = Alma.download_and_extract_files(list(set(m83urls['URL'])))
m83files = m83files

Simbad.add_votable_fields('rvel')
m83simbad = Simbad.query_object('M83')
rvel = m83simbad['RVel_Rvel'][0]*u.Unit(m83simbad['RVel_Rvel'].unit)

for fn in m83files:
  if 'line' in fn:
    cube = SpectralCube.read(fn)
    # Convert frequencies to their rest frequencies
    frange = u.Quantity([cube.spectral_axis.min(),
      cube.spectral_axis.max()]) * (1+rvel/constants.c)

    # Query the top 20 most common species in the frequency range of the
    # cube with an upper energy state <= 50K
    lines = Splatalogue.query_lines(frange[0], frange[1], top20='top20',
      energy_max=50, energy_type='eu_k',
      only_NRAO_recommended=True)

    lines.pprint()

    # Change the cube coordinate system to be in velocity with respect
    # to the rest frequency (in the M83 rest frame)
    rest_frequency = lines['Freq-GHz'][0]*u.GHz / (1+rvel/constants.c)
    vcube = cube.with_spectral_unit(u.km/u.s,
      rest_value=rest_frequency, velocity_convention='radio')

    # Write the cube with the specified line name
    fmt = "{Species}{Resolved QNs}"
    row = lines[0]
    linename = fmt.format(**dict(zip(row.colnames,row.data)))
    vcube.write('M83_ALMA_{linename}.fits'.format(linename=linename))
```

## Example #7

Find ALMA pointings that have been observed toward M83, then overplot
the various fields-of view on a 2MASS image retrieved from SkyView. See
[http://nbviewer.ipython.org/gist/keflavich/19175791176e8d1fb204][25]
for the notebook. There is an even more sophisticated version at
[http://nbviewer.ipython.org/gist/keflavich/bb12b772d6668cf9181a][26],
which shows Orion KL in all observed bands.

```python

# Querying ALMA archive for M83 pointings and plotting
# them on a 2MASS image

In [2]: from astroquery.alma import Alma
        from astroquery.skyview import SkyView
        import string
        from astropy import units as u
        import pylab as pl
        import aplpy

# Retrieve M83 2MASS K-band image:

In[3]: m83_images = SkyView.get_images(position='M83', survey=['2MASS-K'], pixels=1500)

# Retrieve ALMA archive information *including* private data and
# non-science fields:

In[4]: m83 = Alma.query_object('M83', public=False, science=False)

In[5]: m83

# Parse components of the ALMA data.  Specifically, find the frequency
# support - the frequency range covered - and convert that into a central
# frequency for beam radius estimation.

In[6]: def parse_frequency_support(frequency_support_str):
            supports = frequency_support_str.split("U")
            freq_ranges = [(float(sup.strip('[] ').split("..")[0]),
              float(sup.strip('[] ').split("..")[1].split(',')[0].strip(string.letters)))
              *u.Unit(sup.strip('[] ').split("..")[1].split(',')[0].strip(string.punctuation+string.digits)) for sup in supports]
              for sup in supports]
            return u.Quantity(freq_ranges)

       def approximate_primary_beam_sizes(frequency_support_str):
          freq_ranges = parse_frequency_support(frequency_support_str)
          beam_sizes = [(1.22*fr.mean().to(u.m, u.spectral())/(12*u.m)).to(u.arcsec,u.dimensionless_angles())
            for fr in freq_ranges]
          return u.Quantity(beam_sizes)

In[7]: primary_beam_radii = [approximate_primary_beam_sizes(row['Frequency support']) for row in m83]

# Compute primary beam parameters for the public and private
# components of the data for plotting below.

In[8]: print "The bands used include: ",np.unique(m83['Band'])

In[9]: private_circle_parameters = [(row['RA'],row['Dec'],np.mean(rad).to(u.deg).value)
          for row,rad in zip(m83, primary_beam_radii)
          if row['Release date']!='' and row['Band']==3]
       public_circle_parameters = [(row['RA'],row['Dec'],np.mean(rad).to(u.deg).value)
          for row,rad in zip(m83, primary_beam_radii)
          if row['Release date']=='' and row['Band']==3]
       unique_private_circle_parameters = np.array(list(set(private_circle_parameters)))
       unique_public_circle_parameters = np.array(list(set(public_circle_parameters)))

       print "BAND 3"
       print "PUBLIC:  Number of rows: {0}.  Unique pointings: {1}".format(len(m83), len(unique_public_circle_parameters))
       print "PRIVATE: Number of rows: {0}.  Unique pointings: {1}".format(len(m83), len(unique_private_circle_parameters))

       private_circle_parameters_band6 = [(row['RA'],row['Dec'],np.mean(rad).to(u.deg).value)
          for row,rad in zip(m83, primary_beam_radii)
          if row['Release date']!='' and row['Band']==6]
       public_circle_parameters_band6 = [(row['RA'],row['Dec'],np.mean(rad).to(u.deg).value)
          for row,rad in zip(m83, primary_beam_radii)
          if row['Release date']=='' and row['Band']==6]

# Show all of the private observation pointings that have been acquired
In[10]: fig = aplpy.FITSFigure(m83_images[0])
        fig.show_grayscale(stretch='arcsinh')
        fig.show_circles(unique_private_circle_parameters[:,0],
          unique_private_circle_parameters[:,1],
          unique_private_circle_parameters[:,2],
          color='r', alpha=0.2)

# In principle, all of the pointings shown below should be downloadable
# from the archive:

In[11]: fig = aplpy.FITSFigure(m83_images[0])
        fig.show_grayscale(stretch='arcsinh')
        fig.show_circles(unique_public_circle_parameters[:,0],
          unique_public_circle_parameters[:,1],
          unique_public_circle_parameters[:,2],
          color='b', alpha=0.2)

# Use pyregion to write the observed regions to disk.  Pyregion has a
# very awkward API; there is (in principle) work in progress to improve
# that situation but for now one must do all this extra work.

In[16]: import pyregion
        from pyregion.parser_helper import Shape
        prv_regions = pyregion.ShapeList([Shape('circle',[x,y,r]) for x,y,r in private_circle_parameters])
        pub_regions = pyregion.ShapeList([Shape('circle',[x,y,r]) for x,y,r in public_circle_parameters])
        for r,(x,y,c) in zip(prv_regions+pub_regions, np.vstack([private_circle_parameters, public_circle_parameters])):
          r.coord_format = 'fk5'
          r.coord_list = [x,y,c]
          r.attr = ([], {'color': 'green',  'dash': '0 ',  'dashlist': '8 3 ',  'delete': '1 ',  'edit': '1 ',
            'fixed': '0 ',  'font': '"helvetica 10 normal roman"', 'highlite': '1 ',
            'include': '1 ',  'move': '1 ',  'select': '1 ',  'source': '1',  'text': '',
            'width': '1 '})
        prv_regions.write('M83_observed_regions_private_March2015.reg')
        pub_regions.write('M83_observed_regions_public_March2015.reg')

In[17]: from astropy.io import fits

In[18]: prv_mask = fits.PrimaryHDU(prv_regions.get_mask(m83_images[0][0]).astype('int'),
          header=m83_images[0][0].header)
        pub_mask = fits.PrimaryHDU(pub_regions.get_mask(m83_images[0][0]).astype('int'),
          header=m83_images[0][0].header)

In[19]: pub_mask.writeto('public_m83_almaobs_mask.fits', clobber=True)

In[20]: fig = aplpy.FITSFigure(m83_images[0])
        fig.show_grayscale(stretch='arcsinh')
        fig.show_contour(prv_mask, levels=[0.5,1], colors=['r','r'])
        fig.show_contour(pub_mask, levels=[0.5,1], colors=['b','b'])

# ## More advanced ##
#
# Now we create a 'hit mask' showing the relative depth of each
# observed field in each band

In[21]: hit_mask_band3_public = np.zeros_like(m83_images[0][0].data)
        hit_mask_band3_private = np.zeros_like(m83_images[0][0].data)
        hit_mask_band6_public = np.zeros_like(m83_images[0][0].data)
        hit_mask_band6_private = np.zeros_like(m83_images[0][0].data)

        from astropy import wcs
        mywcs = wcs.WCS(m83_images[0][0].header)

In[22]: for row,rad in zip(m83, primary_beam_radii):
          shape = Shape('circle', (row['RA'], row['Dec'],np.mean(rad).to(u.deg).value))
          shape.coord_format = 'fk5'
          shape.coord_list = (row['RA'], row['Dec'],np.mean(rad).to(u.deg).value)
          shape.attr = ([], {'color': 'green',  'dash': '0 ', 'dashlist': '8 3 ',  'delete': '1 ',
            'edit': '1 ', 'fixed': '0 ',  'font': '"helvetica 10 normal roman"',  'highlite': '1 ',
            'include': '1 ',  'move': '1 ',  'select': '1 ',  'source': '1',  'text': '',
            'width': '1 '})

          if row['Release date']=='' and row['Band']==3:
            (xlo,xhi,ylo,yhi),mask = pyregion_subset(shape, hit_mask_band3_private, mywcs)
            hit_mask_band3_private[ylo:yhi,xlo:xhi] += row['Integration']*mask
          elif row['Release date'] and row['Band']==3:
            (xlo,xhi,ylo,yhi),mask = pyregion_subset(shape, hit_mask_band3_public, mywcs)
            hit_mask_band3_public[ylo:yhi,xlo:xhi] += row['Integration']*mask
          elif row['Release date'] and row['Band']==6:
            (xlo,xhi,ylo,yhi),mask = pyregion_subset(shape, hit_mask_band6_public, mywcs)
            hit_mask_band6_public[ylo:yhi,xlo:xhi] += row['Integration']*mask
          elif row['Release date']=='' and row['Band']==6:
            (xlo,xhi,ylo,yhi),mask = pyregion_subset(shape, hit_mask_band6_private, mywcs)
            hit_mask_band6_private[ylo:yhi,xlo:xhi] += row['Integration']*mask

In[23]: fig = aplpy.FITSFigure(m83_images[0])
        fig.show_grayscale(stretch='arcsinh')
        fig.show_contour(fits.PrimaryHDU(data=hit_mask_band3_public, header=m83_images[0][0].header),
          levels=np.logspace(0,5,base=2, num=6), colors=['r']*6)
        fig.show_contour(fits.PrimaryHDU(data=hit_mask_band3_private, header=m83_images[0][0].header),
          levels=np.logspace(0,5,base=2, num=6), colors=['y']*6)
        fig.show_contour(fits.PrimaryHDU(data=hit_mask_band6_public, header=m83_images[0][0].header),
          levels=np.logspace(0,5,base=2, num=6), colors=['c']*6)
        fig.show_contour(fits.PrimaryHDU(data=hit_mask_band6_private, header=m83_images[0][0].header),
          levels=np.logspace(0,5,base=2, num=6), colors=['b']*6)

In[24]: from astropy import wcs
        import pyregion
        from astropy import log

        def pyregion_subset(region, data, mywcs):
          # Return a subset of an image (`data`) given a region.
          shapelist = pyregion.ShapeList([region])
          if shapelist[0].coord_format not in ('physical','image'):
            # pixel_regions = shapelist.as_imagecoord(self.wcs.celestial.to_header())
            # convert the regions to image (pixel) coordinates
            celhdr = mywcs.sub([wcs.WCSSUB_CELESTIAL]).to_header()
            pixel_regions = shapelist.as_imagecoord(celhdr)
          else:
            # For this to work, we'd need to change the reference pixel
            # after cropping.
            # Alternatively, we can just make the full-sized mask...
            # todo....
            raise NotImplementedError("Can't use non-celestial coordinates with regions.")
            pixel_regions = shapelist
          # This is a hack to use mpl to determine the outer bounds of the regions
          # (but it's a legit hack - pyregion needs a major internal refactor
          # before we can approach this any other way, I think -AG)
          mpl_objs = pixel_regions.get_mpl_patches_texts()[0]

          # Find the minimal enclosing box containing all of the regions
          # (this will speed up the mask creation below)
          extent = mpl_objs[0].get_extents()
          xlo, ylo = extent.min
          xhi, yhi = extent.max
          all_extents = [obj.get_extents() for obj in mpl_objs]
          for ext in all_extents:
            xlo = xlo if xlo < ext.min[0] else ext.min[0]
            ylo = ylo if ylo < ext.min[1] else ext.min[1]
            xhi = xhi if xhi > ext.max[0] else ext.max[0]
            yhi = yhi if yhi > ext.max[1] else ext.max[1]

          log.debug("Region boundaries: ")
          log.debug("xlo={xlo}, ylo={ylo}, xhi={xhi},
            yhi={yhi}".format(xlo=xlo, ylo=ylo, xhi=xhi, yhi=yhi))

          subwcs = mywcs[ylo:yhi, xlo:xhi]
          subhdr = subwcs.sub([wcs.WCSSUB_CELESTIAL]).to_header()
          subdata = data[ylo:yhi, xlo:xhi]

          mask = shapelist.get_mask(header=subhdr,
            shape=subdata.shape)
          log.debug("Shapes: data={0}, subdata={2}, mask={1}".format(data.shape, mask.shape, subdata.shape))
          return (xlo,xhi,ylo,yhi),mask
```

## Example #8

Retrieve data from a particular co-I or PI from the ESO archive

```python
from astroquery.eso import Eso

# log in so you can get proprietary data
Eso.login('aginsburg')
# make sure you don't filter out anything
Eso.ROW_LIMIT = 1e6

# List all of your pi/co projects
all_pi_proj = Eso.query_instrument('apex', pi_coi='ginsburg')

# Have a look at the project IDs only
print(set(all_pi_proj['APEX Project ID']))
# set(['E-095.F-9802A-2015', 'E-095.C-0242A-2015', 'E-093.C-0144A-2014'])

# The full project name includes prefix and suffix
full_proj = 'E-095.F-9802A-2015'
proj_id = full_proj[2:-6]

# Then get the APEX quicklook "reduced" data
tbl = Eso.query_apex_quicklooks(prog_id=proj_id)

# and finally, download it
files = Eso.retrieve_data(tbl['Product ID'])

# then move the files to your local directory
# note that there is no .TAR suffix... not sure why this is
import shutil
for fn in files:
  shutil.move(fn+'.TAR','.')
```

# Available Services

The following modules have been completed using a common API:

* [SIMBAD Queries (astroquery.simbad)][30]
* [VizieR Queries (astroquery.vizier)][31]
* [ESASky Queries (astroquery.esasky)][32]
* [IRSA Dust Extinction Service Queries (astroquery.irsa_dust)][33]
* [NED Queries (astroquery.ned)][34]
* [Splatalogue Queries (astroquery.splatalogue)][35]
* [Vamdc Queries (astroquery.vamdc)][36]
* [IRSA Image Server program interface (IBE) Queries (astroquery.ibe)][37]
* [IRSA Queries (astroquery.irsa)][38]
* [UKIDSS Queries (astroquery.ukidss)][39]
* [1MAGPIS Queries (astroquery.magpis)][40]
* [NRAO Queries (astroquery.nrao)][41]
* [Besancon Queries (astroquery.besancon)][42]
* [NIST Queries (astroquery.nist)][43]
* [NVAS Queries (astroquery.nvas)][44]
* [GAMA Queries (astroquery.gama)][45]
* [ESO Queries (astroquery.eso)][46]
* [xMatch Queries (astroquery.xmatch)][47]
* [Atomic Line List (astroquery.atomic)][48]
* [ALMA Queries (astroquery.alma)][49]
* [Skyview Queries (astroquery.skyview)][50]
* [NASA ADS Queries (astroquery.nasa_ads)][51]
* [HEASARC Queries (astroquery.heasarc)][52]

These others are functional, but do not follow a common & consistent
API:

* [Fermi Queries (astroquery.fermi)][53]
* [SDSS Queries (astroquery.sdss)][54]
* [ALFALFA Queries (astroquery.alfalfa)][55]
* [Spitzer Heritage Archive (astroquery.sha)][56]
* [LAMDA Queries (astroquery.lamda)][57]
* [OGLE Queries (astroquery.ogle)][58]
* [Open Exoplanet Catalogue(astroquery.open_exoplanet_catalogue)][59]
* [CosmoSim Queries (astroquery.cosmosim)][60]
* [HITRAN Queries (astroquery.hitran)][61]

# Catalog, Archive, and Other

A second index of the services by the type of data they serve. Some
services perform many tasks and are listed more than once.

## Catalogs

The first serve catalogs, which generally return one row of information
for each source (though they may return many catalogs that each have one
row for each source)

* [ALFALFA Queries (astroquery.alfalfa)][62]
* [GAMA Queries (astroquery.gama)][63]
* [IRSA Image Server program interface (IBE) Queries (astroquery.ibe)][64]
* [IRSA Queries (astroquery.irsa)][65]
* [IRSA Dust Extinction Service Queries (astroquery.irsa_dust)][66]
* [NED Queries (astroquery.ned)][67]
* [OGLE Queries (astroquery.ogle)][68]
* [Open Exoplanet Catalogue(astroquery.open_exoplanet_catalogue)][69]
* [SDSS Queries (astroquery.sdss)][70]
* [Spitzer Heritage Archive (astroquery.sha)][71]
* [SIMBAD Queries (astroquery.simbad)][72]
* [UKIDSS Queries (astroquery.ukidss)][73]
* [VizieR Queries (astroquery.vizier)][74]
* [xMatch Queries (astroquery.xmatch)][75]

## Archives

Archive services provide data, usually in FITS images or spectra. They
will generally return a table listing the available data first.

* [ALFALFA Queries (astroquery.alfalfa)][76]
* [ALMA Queries (astroquery.alma)][77]
* [ESO Queries (astroquery.eso)][78]
* [Fermi Queries (astroquery.fermi)][79]
* [HEASARC Queries (astroquery.heasarc)][80]
* [IRSA Image Server program interface (IBE) Queries (astroquery.ibe)][81]
* [IRSA Queries (astroquery.irsa)][82]
* [MAGPIS Queries (astroquery.magpis)][83]
* [NED Queries (astroquery.ned)][84]
* [NRAO Queries (astroquery.nrao)][85]
* [NVAS Queries (astroquery.nvas)][86]
* [SDSS Queries (astroquery.sdss)][87]
* [Spitzer Heritage Archive (astroquery.sha)][88]
* [UKIDSS Queries (astroquery.ukidss)][89]
* [Skyview Queries (astroquery.skyview)][90]

## Simulations

Simulation services query databases of simulated or synthetic data

* [Besancon Queries (astroquery.besancon)][91]
* [CosmoSim Queries (astroquery.cosmosim)][92]

## Other

There are other astronomically significant services, e.g. line list and
atomic/molecular cross section and collision rate services, that don’t
fit the above categories.

* [Atomic Line List (astroquery.atomic)][93]
* [LAMDA Queries (astroquery.lamda)][94]
* [NIST Queries (astroquery.nist)][95]
* [Splatalogue Queries (astroquery.splatalogue)][96]
* [NASA ADS Queries (astroquery.nasa_ads)][97]
* [Vamdc Queries (astroquery.vamdc)][98]
* [HITRAN Queries (astroquery.hitran)][99]

## Developer documentation

The [Astroquery API Specification][100] is intended to be kept as consistent as
possible, such that any web service can be used with a minimal learning
curve imposed on the user.

* [Astroquery API Specification][100]
* [Template Module][101]
* [Astroquery Testing][102]

The following Astroquery modules are mostly meant for internal use of
services in Astroquery, you can use them for your scripts, but we don’t
guarantee API stability.

* [Astroquery utils (astroquery.utils)][103]
* [Astroquery query (astroquery.query)][104]

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
[24]: http://astroquery.readthedocs.io/en/latest/api/astroquery.simbad.SimbadClass.html#astroquery.simbad.SimbadClass.list_votable_fields
[25]: http://nbviewer.ipython.org/gist/keflavich/19175791176e8d1fb204
[26]: http://nbviewer.ipython.org/gist/keflavich/bb12b772d6668cf9181a
[30]: http://astroquery.readthedocs.io/en/latest/simbad/simbad.html
[31]: http://astroquery.readthedocs.io/en/latest/vizier/vizier.html
[32]: http://astroquery.readthedocs.io/en/latest/esasky/esasky.html
[33]: http://astroquery.readthedocs.io/en/latest/irsa/irsa_dust.html
[34]: http://astroquery.readthedocs.io/en/latest/ned/ned.html
[35]: http://astroquery.readthedocs.io/en/latest/splatalogue/splatalogue.html
[36]: http://astroquery.readthedocs.io/en/latest/vamdc/vamdc.html
[37]: http://astroquery.readthedocs.io/en/latest/ibe/ibe.html
[38]: http://astroquery.readthedocs.io/en/latest/irsa/irsa.html
[39]: http://astroquery.readthedocs.io/en/latest/ukidss/ukidss.html
[40]: http://astroquery.readthedocs.io/en/latest/magpis/magpis.html
[41]: http://astroquery.readthedocs.io/en/latest/nrao/nrao.html
[42]: http://astroquery.readthedocs.io/en/latest/besancon/besancon.html
[43]: http://astroquery.readthedocs.io/en/latest/nist/nist.html
[44]: http://astroquery.readthedocs.io/en/latest/nvas/nvas.html
[45]: http://astroquery.readthedocs.io/en/latest/gama/gama.html
[46]: http://astroquery.readthedocs.io/en/latest/eso/eso.html
[47]: http://astroquery.readthedocs.io/en/latest/xmatch/xmatch.html
[48]: http://astroquery.readthedocs.io/en/latest/atomic/atomic.html
[49]: http://astroquery.readthedocs.io/en/latest/alma/alma.html
[50]: http://astroquery.readthedocs.io/en/latest/skyview/skyview.html
[51]: http://astroquery.readthedocs.io/en/latest/nasa_ads/nasa_ads.html
[52]: http://astroquery.readthedocs.io/en/latest/heasarc/heasarc.html
[53]: http://astroquery.readthedocs.io/en/latest/fermi/fermi.html
[54]: http://astroquery.readthedocs.io/en/latest/sdss/sdss.html
[55]: http://astroquery.readthedocs.io/en/latest/alfalfa/alfalfa.html
[56]: http://astroquery.readthedocs.io/en/latest/sha/sha.html
[57]: http://astroquery.readthedocs.io/en/latest/lamda/lamda.html
[58]: http://astroquery.readthedocs.io/en/latest/ogle/ogle.html
[59]: http://astroquery.readthedocs.io/en/latest/open_exoplanet_catalogue/open_exoplanet_catalogue.html
[60]: http://astroquery.readthedocs.io/en/latest/cosmosim/cosmosim.html
[61]: http://astroquery.readthedocs.io/en/latest/hitran/hitran.html
[62]: http://astroquery.readthedocs.io/en/latest/alfalfa/alfalfa.html
[63]: http://astroquery.readthedocs.io/en/latest/gama/gama.html
[64]: http://astroquery.readthedocs.io/en/latest/ibe/ibe.html
[65]: http://astroquery.readthedocs.io/en/latest/irsa/irsa.html
[66]: http://astroquery.readthedocs.io/en/latest/irsa/irsa_dust.html
[67]: http://astroquery.readthedocs.io/en/latest/irsa/irsa_dust.html
[68]: http://astroquery.readthedocs.io/en/latest/ogle/ogle.html
[69]: http://astroquery.readthedocs.io/en/latest/open_exoplanet_catalogue/open_exoplanet_catalogue.html
[70]: http://astroquery.readthedocs.io/en/latest/sdss/sdss.html
[71]: http://astroquery.readthedocs.io/en/latest/sha/sha.html
[72]: http://astroquery.readthedocs.io/en/latest/simbad/simbad.html
[73]: http://astroquery.readthedocs.io/en/latest/ukidss/ukidss.html
[74]: http://astroquery.readthedocs.io/en/latest/vizier/vizier.html
[75]: http://astroquery.readthedocs.io/en/latest/xmatch/xmatch.html
[76]: http://astroquery.readthedocs.io/en/latest/alfalfa/alfalfa.html
[77]: http://astroquery.readthedocs.io/en/latest/alma/alma.html
[78]: http://astroquery.readthedocs.io/en/latest/eso/eso.html
[79]: http://astroquery.readthedocs.io/en/latest/fermi/fermi.html
[80]: http://astroquery.readthedocs.io/en/latest/heasarc/heasarc.html
[81]: http://astroquery.readthedocs.io/en/latest/ibe/ibe.html
[82]: http://astroquery.readthedocs.io/en/latest/irsa/irsa.html
[83]: http://astroquery.readthedocs.io/en/latest/magpis/magpis.html
[84]: http://astroquery.readthedocs.io/en/latest/ned/ned.html
[85]: http://astroquery.readthedocs.io/en/latest/nrao/nrao.html
[86]: http://astroquery.readthedocs.io/en/latest/nvas/nvas.html
[87]: http://astroquery.readthedocs.io/en/latest/sdss/sdss.html
[88]: http://astroquery.readthedocs.io/en/latest/sha/sha.html
[89]: http://astroquery.readthedocs.io/en/latest/ukidss/ukidss.html
[90]: http://astroquery.readthedocs.io/en/latest/skyview/skyview.html
[91]: http://astroquery.readthedocs.io/en/latest/besancon/besancon.html
[92]: http://astroquery.readthedocs.io/en/latest/cosmosim/cosmosim.html
[93]: http://astroquery.readthedocs.io/en/latest/atomic/atomic.html
[94]: http://astroquery.readthedocs.io/en/latest/lamda/lamda.html
[95]: http://astroquery.readthedocs.io/en/latest/nist/nist.html
[96]: http://astroquery.readthedocs.io/en/latest/splatalogue/splatalogue.html
[97]: http://astroquery.readthedocs.io/en/latest/nasa_ads/nasa_ads.html
[98]: http://astroquery.readthedocs.io/en/latest/vamdc/vamdc.html
[99]: http://astroquery.readthedocs.io/en/latest/hitran/hitran.html
[100]: http://astroquery.readthedocs.io/en/latest/api.html
[101]: http://astroquery.readthedocs.io/en/latest/template.html
[102]: http://astroquery.readthedocs.io/en/latest/testing.html
[103]: http://astroquery.readthedocs.io/en/latest/utils.html
[104]: http://astroquery.readthedocs.io/en/latest/query.html
