This tool called polpredict.cpp to reconstruct harmonics for a tidal prediction for given times. It is demonstrated for a port or for a region of model output. However the algorithm iterates over each target location rather than being able to accept arrays.

I learnt how to invoke C in a python wrapping using the code that Lovro built to generate velocity map tiles from the anyTide server.

Notes for building anyTide_Cwrapper.py 

	module load anaconda/2.1.0 
	conda create --name anyTide_env matplotlib=1.4.2 numpy=1.9.1 pandas simplejson netCDF4
	source activate anyTide_env

	cd /login/jelt/python/ipynb/anytide/pyTide

	python setup.py install

	python anyTide_Cwrapper.py
