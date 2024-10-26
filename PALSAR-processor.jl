### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 2e1baf1a-93ca-11ef-3482-6184fb11d8fc
begin
	using ImageView
	using Serialization
	using FFTW
end

# ╔═╡ 4ed04bf8-eeee-4c2b-87d1-4517c9bfcc6c
begin
	pathname = "Cape/"    # relative path of folder containing L1.0 data
	imagename = "IMG-HH"  # which image polarization to use - IMG-HH or IMG-HV for PALSAR
end

# ╔═╡ dc7fbf03-372d-47b2-9423-44dd535dc39c
md"## Schema"

# ╔═╡ f0ccf9eb-fa4e-4887-829e-bf56d24a49e1
function bytesToString(file, start, stop)
    len = stop-start+1
    seek(file,start-1)
    String(read(file,len))
end

# ╔═╡ f7774d31-8026-4e7d-a5ef-e902b6f0b67a
function bytesToUInt32(file, start, stop)
    #TODO check if range matches size
    seek(file, start-1)
    num = try
        read(file,UInt32)
    catch
        0
    end
    Int32(ntoh(num))
end

# ╔═╡ 9d466ec7-67c4-4041-97f0-ccbb192178a4
function bytesToUInt16(file, start, stop)
    #TODO check if range matches size
    seek(file, start-1)
    Int16(ntoh(read(file,UInt16)))
end

# ╔═╡ 1a7f78b7-0470-4a8a-8ecf-b6415397f7fc
function bytesToArray(file, start, stop)
    len = stop-start+1
    seek(file, start-1)
    read(file,len)
end

#Image File Descriptor

# ╔═╡ a2c9dde3-9779-430f-87f1-264c2872b44c
imageFileDescrictorScheme = [    
        ("length", (9, 12), bytesToUInt32),
        ("numSignals", (181, 186), bytesToString),
        ("sarDataBytes", (187, 192), bytesToString),
        ("bitsPerSample", (217, 220), bytesToString),
        ("samplesPerPixel", (221, 224), bytesToString), #could be 1 or 2 (IQ or A)
        #245-272 - borders and interleaving
        ("physicalRecordsPerLine", (273, 274), bytesToString),
        #281-400 - unexplored
        ("sarDataFormat", (401, 428), bytesToString),
        ("sarDataFormatTypeCode", (429,432), bytesToString)
        
        ]

# ╔═╡ 8fe67eea-1bea-4034-9d36-83bbfc67b202
signalDataRecordScheme = [
        ("recordSequenceNumber", (1,4), bytesToUInt32),
        ("lengthOfRecord", (9,12), bytesToUInt32),
        ("sarImageDatalineNumber", (13, 16), bytesToUInt32),
        ("pixelCount", (25, 28), bytesToUInt32),
        #33-56 - instrument settings
        #57-84 - chirp characteristics
        ("PRF (mHz)", (57,60), bytesToUInt32),
        ("Chirp type", (67,68), bytesToUInt16),    #0 for linear
        
        #93-388 - platform info - important!!!
        ("Valid", (97,100),bytesToUInt32),
        ("Electronic Squint 1",(101,104),bytesToUInt32),
        ("Electronic Squint 2",(109,112),bytesToUInt32),
        ("Mechanical Squint",(113,116),bytesToUInt32),
        ("FirstSampleRange (m)",(117,120),bytesToUInt32),
        ("Platform Altitude (m)",(141,144),bytesToUInt32),
        ("Platform Ground Speed (cm/s)",(145,148),bytesToUInt32),
        ("signalData", (413, 10800), bytesToArray)
        ]

# ╔═╡ 22257b03-5b3b-45d2-beac-e4e760a227e9
datasetSummaryRecordScheme = [
    ("wavelength", (501,516), bytesToString),
    ("chirpType", (519,534), bytesToString),
    ("coeffs",(534,694),bytesToString),
    ("samplingRate",(711,726), bytesToString),
    ("pulseLength",(743,758),bytesToString),      #us
    ("PRF", (935,950), bytesToString),            #mHz
    ("doppler", (1415,1654), bytesToString)
]

# ╔═╡ a69a88e6-b776-4df8-abd7-0e1e51fb5205
function parseFile(file,scheme,offset=0)
    keyMap = Dict([])
    for i=1:length(scheme)
        keyMap[scheme[i][1]] = i  #make dictionary mapping field names to indices
    end
    results = []
    for (name, (start, stop), f) in scheme
        len  = stop-start+1
        push!(results, f(file,start+offset,stop+offset))
    end
    return (key=keyMap, fields=results)
end

# ╔═╡ 85f42019-ef9c-4942-84e0-15165c9982a2
md"""
## File Processing

This stage converts a `.0__A` file downloaded from ASF Vertex into a serialized Julia object in a `.ser` file. Can be skipped if you've already produced a `.ser` file for the SAR image you want to process.
"""

# ╔═╡ cf1ce4e3-89f4-4857-922a-1efa34d76e92
md"""
## Chirp Finding
This stage produces a chirp signal that will deconvolve the IQ samples stored in the L1.0 file read in above. Cannot be skipped!
"""

# ╔═╡ 323c60cf-9ad3-4845-a987-32417de259fa
R0 = 848665     #m  
altitude = 628000 #m, nominal
c = 3e8

vorbital = let 
    a = (6371+628)*1000.0  #m            #https://www.eorc.jaxa.jp/ALOS-2/en/about/overview.htm
    G = 6.67408e-11
    Me= 5.97219e24       #kg
    P = sqrt(4*pi^2/(G*Me)*a^3)
    print("Period: ")
    print(round(P/60,digits=4))
    println(" m")
    2pi*a/P            #These are terrible assumptions!!!!!!!
end
vorbital = 7593
Vr = vorbital  # TODO: remove reference to redundant var Vr
println("Orbital Velocity: ",vorbital," m/s")

wavelength = parse(Float32,rec.fields[rec.key[  "wavelength" ]]) #m 
#antenna length
La = 8.9 #m

#beamwidth
bw = 0.886*wavelength/La

# ╔═╡ b66bbf2b-82e5-4de8-9746-57ef6658c57f
nadirAngle = acos(altitude/R0)
println(nadirAngle*180/pi)
RD = 1/2*c*rangeCells/16000000   # TODO where does 16,000,000 come from? 16 MHz sampling rate??
rangeCellLength = 1/rangeCells * (sqrt((R0+RD)^2-altitude^2)-sqrt(R0^2-altitude^2))

# compare the width and height of a pixel (cell)
# to get an approximate scaling for image formation
println("Azimuth:     ", round(vorbital/PRF, digits=3), " m")
println("Range Cell: ",round(rangeCellLength, digits=3), " m")

# approximate "y/x" scaling needed to render images without stretching
aspectRatio = rangeCellLength/(vorbital/PRF)
println("Aspect Ratio: ",round(aspectRatio,digits=1))

# ╔═╡ fe02d877-0df0-4811-99e9-4a6640c81cf3
md"""
## Image Forming
This stage converts the `...-signal-records.ser` object produced above into an early-stage image. The output of these blocks is a range-compressed complex image in a `.rcc` file.
"""

# ╔═╡ d95eaadd-f6ec-4c21-8488-6079963fb946
# TODO: serialize smallSignals?
#Serialization.serialize(open("$pathname/smallSignalsBigger.ser","w"),smallSignals)

# ╔═╡ 4a80c281-68c0-48dd-ba3e-1031e289e70a
md"""
## Loading Range Corrected Complex Image File
Can start here if a range-compressed complex image is already available from IMG-HH.rcc or similar.
"""

# ╔═╡ 5b131478-a2fd-48ba-91f2-83f551d74e12
md"""
### Range Cell Migration Correction
Going to using the Range Doppler Algorithm. Take range compressed data, fourier transform along each line of constant range, then interpolate by a azimuth-frequency-dependent amount to correct for range cell migration.

Could image at this point, but probably won't look different. 

After that, do azimuth compression as usual. To be efficient, don't apply IFFT to RCM corrected data and instead use that direction in the azimuth convolutions.

Range shift at each frequency is:
$$\Delta R(f_n) = \frac{\lambda^2 R_0 f_n^2}{8 V_r^2}$$
Where $R_0$ is the distance of closest approach, $V_r$ is the effective radar velocity (Cummings and Wong pg. 235)

"""

# ╔═╡ 551a1075-7536-44b3-a87e-4c0eab6c14ff
md"""
### Azimuth Compression!!
need PRF, altitude, initial range, ground speed, wavelength
$$R(s) = \sqrt{R_0^2+s^2v^2} = R_0+\dot{R}_0s+\frac{1}{2}\ddot{R}_0s^2$$
$$C(s) = e^{i \frac{4\pi}{\lambda} R(s)}$$ (4pi comes from: phase shift due to distance is 2pi, but you go there and back so phase shift is doubled)
Sample $C(s)$ at $s=n/PRF$

For zero doppler shift (kinda naive case but fine), $\dot{R}_0 = 0$ and $\ddot{R}_0 = \frac{v^2}{R_0}$

**Note:** this all only works for the first pixels!! Need to correct $R_0$ as we move out from the ground track.

"""

# ╔═╡ b74de37d-d437-4023-accf-796c9589851a
begin
	theta(s,R) = atan(vorbital*s/R)
	
	#one way beam pattern:
	p(a) = sinc(a*La/wavelength)
	w(s,R) = p(theta(s,R))^2
	
	Rc = R0+10000        # TODO where does 10000 come from? I think this is a focusing dist choice
	
	R(s) = Rc - 1/2*vorbital^2/Rc*s^2
	
	C(s) = exp(-4pi*im/wavelength*R(s))*w(s,Rc)
	
	complexAzimuthFFT = let
	    width = 200   # TODO where does this come from?
	    s = 1/PRF*(range(-width, stop = width) |> collect)
	
	    sig = C.(s)/sqrt(width)
	    
	    azimuth = vcat(sig, zeros(Complex{Float32},size(cimg)[1]-length(sig)))
	
	    fft(azimuth)
	end
	
	print("Ready")
end

# ╔═╡ d44a11bd-2d77-4ad6-b60d-d9a463874f26
md"## Images"

# ╔═╡ 24eb33eb-3a64-4418-97be-5ae81a902d16
imshow(rawMagnitude);             # show raw echo image


# ╔═╡ fb598d4b-354d-4245-b7ff-e5c387849f00
imshow(rangeCompressedMagnitude); # show range compressed image (chirp deconvolved)


# ╔═╡ ba323a42-5a62-46a3-aabf-c133dd3b8683
imshow(rccftpre);                 # show FFT of deconvolved image w/ RCM curves


# ╔═╡ fe31f73b-4d40-4118-8985-4acbe07ceff3
imshow(rccftpost);                # show FFT of deconvolved image w/ RCM curves corrected


# ╔═╡ ce4aaa78-de6c-4fc7-9b71-ae57c66b48a7
imshow(azcompmag);                # show final azimuth compressed image


# ╔═╡ fe7e4e54-6115-4780-8d2b-684b5fc64e38
imshow(log.(azcompmag));          # show log-scale final image


# ╔═╡ 42512145-c04e-4545-9b75-4eb7f927c149
azcomp = Serialization.deserialize(open("$pathname/$imagename.slc","r"));


# ╔═╡ 90b23695-0354-4c08-bb4b-7ef4b03c567b
# save range and azimuth compressed file as a "single-look complex"
Serialization.serialize(open("$pathname/$imagename.slc","w"),azcomp)

# ╔═╡ 9fda91dd-6f50-48b1-8b26-2e822562916d
# show a subsection of the image at full resolution
imshow(reverse(abs.(view(cimg,(1:4:10000).+16000,(1:1:1600).+1000)),dims=1));

# ╔═╡ a68b1121-cf14-47cd-894c-68b64b27ad78
begin
	# load serialized version of signalRecords if not loaded from file above
	signalRecords = Serialization.deserialize(open("$pathname/$imagename-signal-records.ser","r"))
	print("Loaded signalRecords")
end

# ╔═╡ 9d59b5bd-e93b-4486-8dea-b733d6e19566
begin
	#now we want to shift each frequency in range space, so make each frequency bin a column for speed
	cimg = cimg'
	shape = size(cimg)
	
	slantRes = 1/2*c/sampleRate
	
	for i = 1:shape[2]
	    n = i
	    #the "highest frequencies" are actually the negative frequencies aliased up!
	    if n>shape[2]/2
	        n = shape[2]-i
	    end
	    fn = (n-1)/shape[2]*PRF    #check this but pretty sure
	    
	    #range migration distance in meters
	    ΔR = wavelength^2*R0*fn^2/(8*Vr^2)
	    cellshift = Integer(round(ΔR/slantRes))
	    
	    #interpolation
	    #NEAREST NEIGHBOR - bad!
	    cimg[:,i] = vcat(cimg[cellshift+1:shape[1],i],zeros(Complex{Float64},cellshift))
	    
	    if n%10000==0
	        print(n)
	        print("  ")
	        println(cellshift)
	    end
	end
	cimg = cimg'
	println("Done")
end

# ╔═╡ 4f9be0f5-05fc-46fd-9ed8-428d281ebaf9
# run an fft on each column of cimg (echos are rows here)
cimg = Complex{Float16}.(fft(Complex{Float32}.(cimg),(1)));

# ╔═╡ dec4f34e-b807-4306-97c7-38389a2c4689
begin
	# show fourier transformed cimg
	# the curves that will be corrected
	# by range cell migration should be visible
	shape = size(cimg)
	rccftpre = (abs.(view(cimg,
	                2:100:shape[1],
	                3600+54:1:3600+473)))
	imshow(rccftpre);
end

# ╔═╡ 5ebaa4a6-4dd6-4c30-b6c1-b22fb34319af
begin
	#sub image formation:
	sampleNum = rangeCells #how many samples of each echo to keep - keep all range cells by default
	echoNum = 35000  #how many echos to keep
	echostart = 1
	
	# smallSignals is matrix containing a subset (defined by sampleNum and echoNum)
	# of the echo signals in signalRecords
	smallSignals = zeros(Complex{Float16},sampleNum,echoNum)
	
	#This is soooooo much faster than hcatS!!!!
	#https://stackoverflow.com/questions/38308053/julia-how-to-fill-a-matrix-row-by-row-in-julia
	for i = echostart+1:echostart+echoNum
	    line = signalRecords[i]
	    I = Float16.(line.fields[line.key["I"]])
	    Q = Float16.(line.fields[line.key["Q"]])
	    
	    smallSignals[:,i-echostart] = Complex.(I[1:sampleNum],Q[1:sampleNum])
	end
	
	signalRecords = []  #free up signalRecords - now using smallSignals
	
	print("Done")
end

# ╔═╡ a4b94eea-6f4e-406c-96de-4b62cd3636fe
begin
	# deconvolution by the chirp signal
	
	shape = size(smallSignals)
	
	# add zero padding at the beginning of each pulse echo (each column is an echo)
	cimg = vcat(zeros(Complex{Float32},(pulseSamples,shape[2])),
	             Complex{Float32}.(smallSignals));
	
	fft!(cimg,(1)); # perform an FFT on each column (each pulse echo)
	
	cimg =  cimg .* conj.(chirpFFT) ; # convolution with chirp signal performed in frequency domain
	
	ifft!(cimg,(1)); # perform an inverse FFT on each column (each pulse echo)
	cimg = Complex{Float16}.(cimg');  #transpose so that each echo is a horizontal line
	smallSignals = [] # free up memory of smallSignals
end

# ╔═╡ f1808c4f-fd2e-45cc-9f3f-5eb387fa9cfd
begin
	file = open("$pathname/LED.0__A")
	rec = parseFile(file,datasetSummaryRecordScheme,720)
	close(file)
	
	PRF = parse(Float64,rec.fields[rec.key["PRF"]])/1000
	rangeCells = 5000
	
	
	sampleRate = parse(Float64,rec.fields[rec.key["samplingRate"]])*1e6     #Hz
	pulseSamples = let
	    pulseLength = parse(Float64,rec.fields[rec.key["pulseLength"]])*1e-6   #s
	    #PRF = parse(Float64,rec.fields[rec.key["PRF"]])/1000
	    Integer(floor(pulseLength*sampleRate))
	end
	
	#process header
	chirpFFT = let
	    f = parse(Float64,split(rec.fields[rec.key["coeffs"]])[1])
	    fdot = -parse(Float64,split(rec.fields[rec.key["coeffs"]])[2])
	    sampleRate = parse(Float64,rec.fields[rec.key["samplingRate"]])*1e6     #Hz
	    pulseLength = parse(Float64,rec.fields[rec.key["pulseLength"]])*1e-6   #s
	    #PRF = parse(Float64,rec.fields[rec.key["PRF"]])/1000
	    pulseSamples = Integer(floor(pulseLength*sampleRate))
	    #pulseSamples = 400
	    print("Max IF freq: ")
	    println(fdot*pulseLength)
	
	    Sif(t) = exp(pi*im*fdot*t^2)
	
	    t = 1/sampleRate* (range(1, stop = pulseSamples) |> collect)
	
	    t = t.-maximum(t)/2
	
	    sig = Sif.(t)/sqrt(pulseSamples)
	
	    chirp = vcat(sig,zeros(Complex{Float32}, rangeCells))
	
	    fft(chirp)
	end
	
	print("Ready")
end

# ╔═╡ 0604faa5-c48a-4b5b-84c2-15f05785238f
begin
	# show fourier transformed cimg
	# the curves that will be corrected
	# by range cell migration should be visible
	shape = size(cimg)
	rccftpost = (abs.(view(cimg,
	                2:100:shape[1],
	                3600+54:1:3600+473)))
	imshow(rccftpost);
end

# ╔═╡ d1968e2c-e012-428b-bb19-df4b4ea7c311
begin
	for i = 1:size(cimg)[2]
	    line = Complex{Float32}.(cimg[:,i])
	    
	    ####### Azimuth Compression
	    
	    #lineFFT = fft(line)
	    lineFFT = cimg[:,i]
	    #lineFFT = fft(Complex{Float32}.(cimg[:,i]))
	    crossCorrelated = AbstractFFTs.ifft(conj.(complexAzimuthFFT).*lineFFT)
	    
	    ####### End Azimuth Compression
	    
	    
	    #complex = Complex.(I,Q)
	    result = abs.( crossCorrelated )
	    cimg[:,i] = result
	    if i%1000 == 1
	        print("#")
	    end
	end
	
	shape = size(cimg)
	azcompmag = abs.(view(cimg,1:16:shape[1],1:4:shape[2]));
	azcompmag = reverse(azcompmag,dims=1)
	imshow(azcompmag);
end

# ╔═╡ c855ab71-65e6-4805-8d09-ba8b028f80af
begin
	cimg = Serialization.deserialize(open("$pathname/$imagename.rcc","r"))
	print("Loaded")
end

# ╔═╡ 5780e24c-8cf2-403b-af60-8f1f8134cd20
begin
	shape = size(smallSignals)
	rawMagnitude = abs.(view(smallSignals,1:10:shape[1],1:40:shape[2]));
	rawMagnitude = reverse(rawMagnitude,dims=1)
	imshow(rawMagnitude);
	# TODO: might want to extend the width of the image by pulseSamples before downsizing,
	# as was done in the original version. See snippet below:
	# vcat(zeros(Complex{Float16},pulseSamples),smallSignals[:,i])
end

# ╔═╡ 8c0bd2ab-8b14-4375-bf7f-42e1e94dd6e4
begin
	file = open("$pathname/$imagename.0__A")
	
	imageDescriptor = parseFile(file,imageFileDescrictorScheme)
	sarDataBytes    = parse(Int64,imageDescriptor.fields[imageDescriptor.key["sarDataBytes"]])
	numSignals      = parse(Int64,imageDescriptor.fields[imageDescriptor.key["numSignals"]])
	
	signalRecords = [ ]
	
	for i = 1:(numSignals-1)
	    parsedRecord = parseFile(file,signalDataRecordScheme,720+sarDataBytes*i)
	    signal = parsedRecord.fields[parsedRecord.key["signalData"]]
	    
	    #TODO  - is 15 or 16 better?? compute mean of signals?
	    #TODO  - DEAL WITH VALUES > 0x1f!!
	    signal = map(x->min(0x1f,x), signal)
	    
	    signal = Int8.(signal)-Int8.(15*ones(size(signal)))
	    I = signal[1:2:length(signal)]
	    Q = signal[2:2:length(signal)]
	
	    parsedRecord.key["I"] = length(parsedRecord.fields)+1;
	    parsedRecord.key["Q"] = length(parsedRecord.fields)+2;
	    push!(parsedRecord.fields, I)
	    push!(parsedRecord.fields, Q)
	    
	    parsedRecord.fields[parsedRecord.key["signalData"]] = []  #remove original signalData from record!
	    
	    push!(signalRecords,parsedRecord)
	    if(i%10000==0)
	        print(i)
	        print("   ")
	        print(Base.summarysize(parsedRecord))
	        print("B    ")
	        print(Base.summarysize(signalRecords)/1e6)
	        println(" MB")
	    end
	end
	print(Base.summarysize(signalRecords)/1e6)
	println(" MB")
	
	close(file)
	
	# save signalRecords for later use
	Serialization.serialize(open("$pathname/$imagename-signal-records.ser","w"),signalRecords)
	signalRecords = [];  # free up the memory of the signalRecords object - it'll be read back in later
end

# ╔═╡ bf110ea9-0741-4d7a-aaa8-5f582f7b7b25
begin
	Serialization.serialize(open("$pathname/$imagename.rcc","w"),cimg)
	cimg = [];
	#GC.gc();
end

# ╔═╡ 85355608-83af-4f20-8bdf-4d485aed6578
begin
	shape = size(cimg)
	rangeCompressedMagnitude = abs.(view(cimg,1:40:shape[1],1:10:shape[2]));
	rangeCompressedMagnitude = reverse(rangeCompressedMagnitude,dims=1)
	imshow(rangeCompressedMagnitude);
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
ImageView = "86fae568-95e7-573e-a6b2-d8a6b900c9ef"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[compat]
FFTW = "~1.8.0"
ImageView = "~0.12.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "90be92a6339a834650603c2111f0dd41f24aee04"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+2"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "7b6ad8c35f4bc3bca8eb78127c8b99719506a5fb"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GTK4_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "Graphene_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Libepoxy_jll", "Libtiff_jll", "Pango_jll", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libX11_jll", "Xorg_libXcursor_jll", "Xorg_libXdamage_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "Xorg_libXrender_jll", "gdk_pixbuf_jll", "iso_codes_jll", "xkbcommon_jll"]
git-tree-sha1 = "b13d8bc60311f88b887f7886633fc396134af050"
uuid = "6ebb71f1-8434-552f-b6b1-dc18babcca63"
version = "4.14.5+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphene_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "61850a17f562453e3485a489c9c8cccb3abcab93"
uuid = "75302f13-0b7e-5bab-a6d1-23fa92e4c2ea"
version = "1.10.6+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Gtk4]]
deps = ["BitFlags", "CEnum", "Cairo", "Cairo_jll", "ColorTypes", "FixedPointNumbers", "GTK4_jll", "Glib_jll", "Graphene_jll", "Graphics", "JLLWrappers", "Libdl", "Librsvg_jll", "Pango_jll", "Preferences", "Reexport", "Scratch", "Xorg_xkeyboard_config_jll", "adwaita_icon_theme_jll", "gdk_pixbuf_jll", "hicolor_icon_theme_jll", "libpng_jll"]
git-tree-sha1 = "f8500791b0cc02587aaa1456b1d8c17c42d030ab"
uuid = "9db2cae5-386f-4011-9d63-a5602296539b"
version = "0.7.0"

[[deps.GtkObservables]]
deps = ["Cairo", "Colors", "Dates", "FixedPointNumbers", "Graphics", "Gtk4", "IntervalSets", "LinearAlgebra", "Observables", "PrecompileTools", "Reexport", "RoundingIntegers"]
git-tree-sha1 = "001628258ac5908ea87e754f4bfb9fc318ddbaff"
uuid = "8710efd8-4ad6-11eb-33ea-2d5ceb25a41c"
version = "2.1.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "b2a7eaa169c13f5bcae8131a83bc30eff8f71be0"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.2"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.ImageView]]
deps = ["AxisArrays", "Cairo", "Compat", "Graphics", "Gtk4", "GtkObservables", "ImageBase", "ImageCore", "ImageMetadata", "MultiChannelColors", "PrecompileTools", "RoundingIntegers", "StatsBase"]
git-tree-sha1 = "e535c709a4fb6f0ae65026fa9bfac624abec016f"
uuid = "86fae568-95e7-573e-a6b2-d8a6b900c9ef"
version = "0.12.6"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libepoxy_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "7a0158b71f8be5c771e7a273183b2d0ac35278c5"
uuid = "42c93a91-0102-5b3f-8f9d-e41de60ac950"
version = "1.5.10+0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultiChannelColors]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "Compat", "FixedPointNumbers", "LinearAlgebra", "Reexport", "Requires"]
git-tree-sha1 = "d29b08ad606124069ca263f9439c00c23f531ea6"
uuid = "d4071afc-4203-49ee-90bc-13ebeb18d604"
version = "0.1.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RoundingIntegers]]
git-tree-sha1 = "99acd97f396ea71a5be06ba6de5c9defe188a778"
uuid = "d5f540fe-1c90-5db3-b776-2e2f362d9394"
version = "1.1.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdamage_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll"]
git-tree-sha1 = "fe4ffb2024ba3eddc862c6e1d70e2b070cd1c2bf"
uuid = "0aeada51-83db-5f97-b67e-184615cfc6f6"
version = "1.1.5+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.adwaita_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "hicolor_icon_theme_jll"]
git-tree-sha1 = "28401767f30e5743ef5e3b0be71417bc911d3952"
uuid = "b437f822-2cd6-5e08-a15c-8bac984d38ee"
version = "43.0.1+0"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "86e7731be08b12fa5e741f719603ae740e16b666"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.10+0"

[[deps.hicolor_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b458a6f6fc2b1a8ca74ed63852e4eaf43fb9f5ea"
uuid = "059c91fe-1bad-52ad-bddd-f7b78713c282"
version = "0.17.0+3"

[[deps.iso_codes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d837a5d2a19d54243dafb6cc98d0f590a603dfa1"
uuid = "bf975903-5238-5d20-8243-bc370bc1e7e5"
version = "4.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═2e1baf1a-93ca-11ef-3482-6184fb11d8fc
# ╠═4ed04bf8-eeee-4c2b-87d1-4517c9bfcc6c
# ╠═dc7fbf03-372d-47b2-9423-44dd535dc39c
# ╠═f0ccf9eb-fa4e-4887-829e-bf56d24a49e1
# ╠═f7774d31-8026-4e7d-a5ef-e902b6f0b67a
# ╠═9d466ec7-67c4-4041-97f0-ccbb192178a4
# ╠═1a7f78b7-0470-4a8a-8ecf-b6415397f7fc
# ╠═a2c9dde3-9779-430f-87f1-264c2872b44c
# ╠═8fe67eea-1bea-4034-9d36-83bbfc67b202
# ╠═22257b03-5b3b-45d2-beac-e4e760a227e9
# ╠═a69a88e6-b776-4df8-abd7-0e1e51fb5205
# ╠═85f42019-ef9c-4942-84e0-15165c9982a2
# ╠═8c0bd2ab-8b14-4375-bf7f-42e1e94dd6e4
# ╠═cf1ce4e3-89f4-4857-922a-1efa34d76e92
# ╠═f1808c4f-fd2e-45cc-9f3f-5eb387fa9cfd
# ╠═323c60cf-9ad3-4845-a987-32417de259fa
# ╠═b66bbf2b-82e5-4de8-9746-57ef6658c57f
# ╠═fe02d877-0df0-4811-99e9-4a6640c81cf3
# ╠═a68b1121-cf14-47cd-894c-68b64b27ad78
# ╠═5ebaa4a6-4dd6-4c30-b6c1-b22fb34319af
# ╠═d95eaadd-f6ec-4c21-8488-6079963fb946
# ╠═5780e24c-8cf2-403b-af60-8f1f8134cd20
# ╠═a4b94eea-6f4e-406c-96de-4b62cd3636fe
# ╠═85355608-83af-4f20-8bdf-4d485aed6578
# ╠═bf110ea9-0741-4d7a-aaa8-5f582f7b7b25
# ╠═4a80c281-68c0-48dd-ba3e-1031e289e70a
# ╠═c855ab71-65e6-4805-8d09-ba8b028f80af
# ╟─5b131478-a2fd-48ba-91f2-83f551d74e12
# ╠═4f9be0f5-05fc-46fd-9ed8-428d281ebaf9
# ╠═dec4f34e-b807-4306-97c7-38389a2c4689
# ╠═9d59b5bd-e93b-4486-8dea-b733d6e19566
# ╠═0604faa5-c48a-4b5b-84c2-15f05785238f
# ╠═551a1075-7536-44b3-a87e-4c0eab6c14ff
# ╠═b74de37d-d437-4023-accf-796c9589851a
# ╠═d1968e2c-e012-428b-bb19-df4b4ea7c311
# ╠═90b23695-0354-4c08-bb4b-7ef4b03c567b
# ╠═d44a11bd-2d77-4ad6-b60d-d9a463874f26
# ╠═24eb33eb-3a64-4418-97be-5ae81a902d16
# ╠═fb598d4b-354d-4245-b7ff-e5c387849f00
# ╠═ba323a42-5a62-46a3-aabf-c133dd3b8683
# ╠═fe31f73b-4d40-4118-8985-4acbe07ceff3
# ╠═ce4aaa78-de6c-4fc7-9b71-ae57c66b48a7
# ╠═fe7e4e54-6115-4780-8d2b-684b5fc64e38
# ╠═42512145-c04e-4545-9b75-4eb7f927c149
# ╠═9fda91dd-6f50-48b1-8b26-2e822562916d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
