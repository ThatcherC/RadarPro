# PALSAR Synthetic Aperature RADAR Processor

**TODO**
- [ ] add download instructions for PALSAR imagery
- [ ] add list of SAR jargon terms
- [ ] use just one big matrix, do operations in place, and save off small images at each stage instead of carrying around giant complex matrices from each stage of processing (`cimg`, `shiftedft`, ...)

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
