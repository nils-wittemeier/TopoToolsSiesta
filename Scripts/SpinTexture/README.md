# 1 Calculating the spin texture (SpinTexture.py)

## 1.1 Required files

Spin texture can be calculated from siesta outputs if either the `HSX`, `TSHS` or `SystemLabel.nc` file was generated. When using the `HSX` file `*.ion`  files are required as well for information about the basis set. 

## 1.2 Usage
usage: SpinTexture.py [-h] infile

Calculate the spin texture of a system.

positional arguments:
  infile      input file. json formatted file with calculation parameters

options:
  -h, --help  show this help message and exit

## 1.3 Input file:
Input file should be written in json format specifying the following fields.

 Field            | Type       | Description 
------------------|------------|-------------
 fdf*             | `string`   | Input file of SIESTA calculation
 out*             | `string`   | Name of output file 
 kpath*           | `dict`     | Specification of reciprocal space sampling (see section 1.3.1)

(*) required field

### 1.3.1 Reciprocal space sampling

 Field            | Type       | Description 
------------------|------------|-------------
 mode             | `string`   | Type of reciprocal space sampling either: [`BandStructure`](https://sisl.readthedocs.io/en/latest/api/generated/sisl.physics.BandStructure.html), [`MonkhorstPack`](https://sisl.readthedocs.io/en/latest/api/generated/sisl.physics.MonkhorstPack.html), [`BrillouinZone`](https://sisl.readthedocs.io/en/latest/api/generated/sisl.physics.BrillouinZone.html) [`Circle`](https://sisl.readthedocs.io/en/latest/api/generated/sisl.physics.BrillouinZone.html#sisl.physics.BrillouinZone.param_circle), `ConstantEnergy`
 kwargs           | `dict`    | Parameters for the sampling method (see below). For a complete list refer to the [sisl documentation](https://zerothi.github.io/sisl/index.html))

(*) required field

 #### BandStructure
 
 Field            | Type                  | Description 
------------------|-----------------------|-------------
 points*          | `float array`         | A list of points that are the corners of the path
 divisions*       | `int` or `int array`  | If integer: total number of points. If array_like: number of points in each segment.
 names            | `str array`           | The associated names of the points on the Brillouin Zone path

#### MonkhorstPack keyword arguments
 
 Field            | Type                     | Description 
------------------|--------------------------|-------------
nkpt*             | `int array`              | A list of number of k-points along each cell direction
displacement      | `float` or `float array` | The displacement of the evenly spaced grid, a single floating number is the displacement for the 3 directions, else they are the individual displacements (Default: None)
centered          | `bool`                   | Whether the k-points are Γ-centered for zero displacement (Default: True)

(*) required field

 #### BrillouinZone keyword arguments

 Field            | Type       | Description 
------------------|------------|-------------
 k                | `list` of `list` of `float`  | List of k points in reciprocal coordinates

(*) required field


#### Circle keyword arguments

 Field            | Type             | Description 
------------------|------------------|-------------
 N_or_dk*         | `int` or `float` | Number of k-points generated using the parameterization (if an integer), otherwise it specifies the discretization length on the circle (in 1/Ang), If the latter case will use less than 4 points a warning will be raised and the number of points increased to 4.
 kR*              | `float`          | Radius of the k-point. In 1/Ang
 normal*          | `float array`    | Normal vector to determine the circle plane
 origin           | `float array`    | Origin of the circle used to generate the circular parameterization

(*) required field

#### Constant Energy

 Field            | Type             | Description 
------------------|------------------|-------------
 energy*          | `float`          | Energy at which the path is supposed to be
 kRmin*           | `float`          | Lower bound for path optimization
 kRmax*           | `float`          | Upper bound for path optimization
 ktol             | `float`          | Tolerance of optimization
 kwrgs            | `dict`           | See keywords arguments for "Circle"

(*) required field

##	1.4 Example: 

### 1.4.1 Spin texture along the M - Γ - K - Γ in a 2D material in a hexagonal material
    {
        "fdf" : "in.fdf" ,
        "out" : "Label.spin.bands.",
        "kpath" : { 
            "mode" : "BandStructure",
            "kwargs":  {
                "points" : [
                    [0.5    , 0.     , 0.],
                    [0.     , 0.     , 0.],
                    [0.33333, 0.33333, 0.],
                    [0.     , 0.     , 0.]
                ],
                "divisions": 300,
                "names" :  ["M", "$\\Gamma$", "K", "M"]
            }
        }
    }

### 1.4.2 Spin texture along a circle in a 2D material
{
    "fdf" : "in.fdf" ,
    "out" : "Label.spin.circle.dat",
    "kpath" : { 
        "mode" : "Circle",
        "kwargs":  {
            "kR" : 0.0085,
            "origin": [0.0    , 0.0    , 0.0],
            "normal": [0.0    , 0.0    , 1.0],
            "N_or_dk": 25,
            "loop": true
        }
    }
}

### 1.4.3 Spin texture along a constant energy path
{
    "fdf" : "in.fdf" ,
    "out" : "Label.spin.circle.dat",
    "kpath" : { 
        "mode" : "constE",
        "energy" : 0.33,
        "etol": 0.01,
        "kRmin" : 0.080,
        "kRmax" : 0.140,
        "kwargs": {
            "origin": [0.0    , 0.0    , 0.0],
            "normal": [0.0    , 0.0    , 1.0],
            "N_or_dk": 49,
            "loop": true
        }
    }
}

# 2. Plot spin texture (PlotSpinTexture.py)

## 2.1 Required files
 - output file generated by SpinTexture.py
 - an input file with plot parameters

## 2.2 Usage
usage: PlotSpinTexture.py [-h] infile

Calculate the spin texture of a system.

positional arguments:
  infile      input file. json formatted file with calculation parameters

options:
  -h, --help  show this help message and exit

## 2.3 Input file

Input file should be written in json format specifying the following fields.

 Field            | Type       | Description 
------------------|------------|-------------
 infile*          | `string`   | Input file of SIESTA calculation
 outfile*         | `string`   | Name of output file. If empty plot is displayed instead using `matplotlibs` GUI of save to a file.
 mode             | `string`   | Switch the plotting mode. Either "bands" for a band structure like plot with spin moments as a color scale, or "quiver" for a scatter plot with spin moments as arrows. Default is `"bands"`.
 grid             | `bool`     | Wether to plot a gird. Default is `False`.
 xrange           | `tuple` of `float` | X-axis range. If `None` `matplotlib` will automatically choose a range. Default is `None`.
 yrange           | `tuple` of `float` | Y-axis range. If `None` `matplotlib` will automatically choose a range. Default is `None`.
 xlabel           | `string` | X-axis label. Default is `""`
 ylabel           | `string` | Y-axis label. Default is `""`

(*) required field

## 2.3 Mode-specific input flags

### 2.3.1 Bands mode

 Field            | Type       | Description 
------------------|------------|-------------
 no-spin          | `bool`     | Only plot bands; Ignore spin moments
 only             | `string` or `None`  | Slect which component of the spin moments to plot. Either "x", "y" or "z" or `None`. If `None` all three are displayed. Default is `None`.

### 2.3.1 Quiver mode

 Field            | Type       | Description 
------------------|------------|-------------
 subplot          | `bool`     | Whether to create subplot for every band. Default is `True`.
 erange           | `tuple` of `float` | Select an energy range. Only band with an average energy in the energy range are display. By default all bands are plotted.
 normal           | `list` of `float`  | Normal vector use to select a identify a plane in k space. Default is `[0,0,1]`.
 line             | `bool`     | Whether a line should be drawn through all k points. Default is 'False'.

## 2.3 Example
    python3 PlotSpinTexture.py b-hex.spin.bands \
      --yrange="-2:2" \
      -o SpinTexture.svg
