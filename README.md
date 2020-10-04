# PALSAR Synthetic Aperature RADAR Processor
In 2018 I wrote a synthetic aperture radar processing program as the final project for a class I was taking called [12.421 - Physical Principles of Remote Sensing](http://student.mit.edu/catalog/m12a.html#12.421). 
With an abundance of spare time in 2020, I have re-worked version of that program and made it ready to share! The code here performs the task of *focusing* synthetic aperture radar (SAR) images - it converts the raw voltage readings from the radio receiver aboard JAXA's ALOS/PALSAR satellite to a (fairly) clear and high-resolution image of the Earth from space. 

Right now the code all works well, but needs more commenting to explain the SAR focusing process. However, I've included my 12.421 [final paper](https://github.com/ThatcherC/RadarPro/blob/master/FinalReport/simple-zipped/SARPaper.pdf) in this repo which does a better job of explaining the details of the process. 

**TODO**
- [x] add download instructions for PALSAR imagery
- [ ] add list of SAR jargon terms
- [x] use just one big matrix, do operations in place, and save off small images at each stage instead of carrying around giant complex matrices from each stage of processing (`cimg`, `shiftedft`, ...)
- [ ] Replace/augment the IPython notebook with a Pluto.jl notebook

### Dependencies
- **Julia:** You'll need a recent version of Julia, available for download [here](https://julialang.org/downloads/) (I'm using Julia 1.5)
- **Julia packages:** The necessary Julia packages can be installed by 
  1. Starting a Julia shell by running `julia` in a terminal
  2. Pressing the `]` key to enter the Julia `pkg` shell
  3. Entering the following line:
  
      `add IJulia, ImageView, Serialization, FFTW`
  4. Done! You can close the Julia shell

- **RADAR samples:** I have included one set of radar recordings here, but more are available for free from ASF Vertex, a NASA website (signup required). See the *Acquiring RADAR Samples* section for details. To use the recording bundled with this repo,
download the file `Cape.zip` from the Assets section of the [Releases](https://github.com/ThatcherC/RadarPro/releases) page.
- **Computer:** This code is a big memory-hungry because of the size of these SAR images. It runs well enough on my 2015 MacBook Pro with 8 GB memory, but only if I close most programs and browser tabs. If you run into Julia kernel crashes during processing, try freeing up some RAM by closing any programs or tabs you aren't using at the moment. The code as-written also saves off a couple intermediate files about 750 MB in size, so make sure you have a couple spare gigabytes of storage space too.
  
### Producing Images
**TODO**
  
run `using IJulia; notebook();` to launch the Jupyter notebook environment and you should

### Acquiring RADAR Samples
This project performs the entire synthetic aperture radar image-forming process
starting from the lowest-level data available: the (downconverted) radio signals
received by the SAR satellite after its pulses bounce of the Earth's surface. This
level of processing is often called Level 1.0 or L1.0, with more processed data products
having names like L1.1, L1.2, and L2.0. 

I got my L1.0 date from the [ASF Vertex Portal](https://vertex.daac.asf.alaska.edu/),
which is run by the [Alaska Satellite Facility](https://asf.alaska.edu/) - this is the
main repository for SAR data of all levels as far as I can tell. 

It's free to download SAR products from the ASF Vertex, but you do need to create an 
account in order to download anything.

I have included one L1.0 product from ALOS PALSAR in this repository
courtesy of JAXA and the Japanese Ministry of Economy, Trade, and Industry (METI), which is permitted
under the [Alaska Satellite Facility/ALOS PALSAR EULA](https://asf.alaska.edu/uncategorized/eula/) 
provided the data is used for peaceful purposes and the attributing is given to JAXA/METI.
The L1.0 product I have done most of my testing on is available [here](https://github.com/ThatcherC/RadarPro/releases/download/initial/Cape.zip) - just unzip it into this repo's 
directory on your computer.
