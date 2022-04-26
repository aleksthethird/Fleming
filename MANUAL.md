# Fleming pipeline users manual

## Quick start

- Get the JGT images from wherever.
- Copy the relevant flats into each night's folder.
- In `main.py` change `raw_base_dir` to point to where you downloaded the JGT fields.
- Change `img_base_dir` to where you want to store the reduced images.
- Run with `python main.py`
- For all 18 fields (as of DR2), will take >6 hours.


## Advanced use

There are a lot of debug options for customisation of the pipeline.
I'll outline a few of the most relevant ones here, but all should have at least some 
documentation in the actual code.
(TODO: run `pydoc` on the codebase).


### Re-running fields

Fields need to be fully processed once.
After this, the light curve generation and adjustment can be skipped and you can choose 
to only run the variable selection.

To do this:
```
Pipeline.run_analysis(config)
```

There are optional parameters, as outlined in [`Pipeline.py`](../pipeline/Pipeline.py).

Simiarly, if you only want to create light curve files, you can run
```
Pipeline.images_to_light_curves(config, start_time)
```

with its own set of optional parameters.

If you want even more control over which parts to run, use the "pipeline breakdown" section
in `interactive.ipynb`.


### Customising the configs

The `Config` objects are designed to be extremely customisable and convenient.
See [`Config.py`](../pipeline/Config.py) for all of the options which can be configured.

You do not need to change the values in `Config.py` unless you want to permanently change the defaults.
To change the values for a single field, you can create a config object like so:

```
config = Config(
	raw_image_dir = "C:/jgtdata",
	n_sets = 5,
	set_size = 100,
	amplitude_score_threshold = 0.5,
)
```

Proof-of-concept of a drop-in use of this pipeline with the remote Napier 3 can be seen in the
"debug" section of [`interactive.ipynb`](../interactive.ipynb).

Out of the box, the configs are suitable for running a field taken with Nebulosity.


### Storing the raw images

The storage location of the images is very flexible as well.
Each `Config` object can be configured to look in one place for the raw images, write reduced
images in another place and then write the working files to yet another place.
By default, these three places are all in the same directory as `main.py`.

This is useful if you have the raw images on a separate drive to the one you're working in.

- `raw_image_dir` - Absolute directory of the raw JGT images to process.
- `image_dir` - Absolute directory of where to store the reduced images.
- `workspace_dir` - Absolute directory of where to store all of the working files and light curves etc.


