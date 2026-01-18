### Listing All Available WhiteboxTools Tools (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

Provides an example of how to list all tools currently available within the WhiteboxTools library, useful for discovering functionality.

```Python
print(wbt.list_tools())
```

--------------------------------

### Running a WhiteboxTools Command with Parameters (Bash)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This example demonstrates how to run a specific WhiteboxTools tool (lidar_info) from the command line. It includes setting the working directory, specifying an input file, and enabling verbose output and GeoKeys processing.

```Bash
./whitebox-tools -r=lidar_info --cd="/path/to/data/" -i=input.las --vlr --geokeys
```

--------------------------------

### Running WhiteboxTools from Command Line

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This example demonstrates how to execute a WhiteboxTools command directly from the terminal. It sets the working directory, specifies the tool to run (`DevFromMeanElev`), and provides input and output file paths. The `-v` flag enables verbose output.

```Shell
./whitebox_tools --wd='/Users/johnlindsay/Documents/data/' --run=DevFromMeanElev --input='DEM clipped.dep' --output='DEV raster.dep' -v
```

--------------------------------

### Initializing WhiteboxTools in Python

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This Python script initializes the WhiteboxTools object, demonstrating how to set the executable path and the working directory. It also shows how to optionally configure verbose mode for tool output. This setup is crucial for interacting with WhiteboxTools from Python applications.

```Python
import os
import sys
from whitebox_tools import WhiteboxTools

wbt = WhiteboxTools()

# If the WhiteboxTools executable file (whitbox_tools.exe) is not in the same
# directory as this script, its path will need to be set, e.g.:
wbt.set_whitebox_dir(os.path.dirname(
    os.path.abspath(__file__)) + "/target/release/")  # or simply wbt.exe_path = ...

# Set the working directory. This is the path to the folder containing the data,
# i.e. files sent to tools as input/output parameters. You don't need to set
# the working directory if you specify full path names as tool parameters.
wbt.work_dir = os.path.dirname(os.path.abspath(__file__)) + "/testdata/"

# Sets verbose mode (True or False). Most tools will suppress output (e.g. updating
# progress) when verbose mode is False. The default is True
# wbt.set_verbose_mode(False) # or simply, wbt.verbose = False
```

--------------------------------

### Printing Specific WhiteboxTools Tool Help (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

Illustrates how to retrieve and print detailed help documentation for a specific WhiteboxTools tool, such as 'ElevPercentile', including its input arguments and example usage. Both CamelCase and snake_case tool names are supported.

```Python
print(wbt.tool_help("ElevPercentile"))
```

--------------------------------

### Running LiDAR Tophat Transform Tool (Shell)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command-line example demonstrates how to execute the `LidarTophatTransform` tool within WhiteboxTools. It specifies the input and output LiDAR files, which are provided as zipped LAS archives, and sets a radius parameter for the transformation. The `--wd` flag sets the working directory, and `-v` enables verbose output.

```Shell
>>./whitebox_tools -r=LidarTophatTransform -v --wd="/path/to/data/"
-i="input.las.zip" -o="output.las.zip" --radius=10.0
```

--------------------------------

### Cloning WhiteboxTools Repository using Git

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command clones the WhiteboxTools source code repository from GitHub to the local machine. It is an alternative to manually downloading and decompressing the zipped source code, allowing users with Git installed to quickly obtain the project files.

```Shell
git clone https://github.com/jblindsay/whitebox-tools.git
```

--------------------------------

### Printing WhiteboxTools General Help (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

Illustrates how to retrieve and print the general help documentation for WhiteboxTools, which provides a listing of all available commands and their basic usage.

```Python
print(wbt.help())
```

--------------------------------

### Running a WhiteboxTools Tool (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

Demonstrates how to execute a WhiteboxTools function, such as `elev_percentile`, using its associated method in Python. It shows how to pass input and output file paths along with tool-specific parameters. An optional custom callback can be provided for processing tool output.

```Python
wbt.elev_percentile("DEM.tif", "output.tif", 15, 15)
```

--------------------------------

### Displaying WhiteboxTools Command-Line Help (Bash)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This command displays the general help information for the WhiteboxTools command-line interface, listing all available flags and their purposes. It is typically run from the terminal after navigating to the WhiteboxTools directory.

```Bash
./whitebox_tools --help
```

--------------------------------

### Compiling WhiteboxTools with Cargo (Release Build)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command uses Cargo, Rust's package manager and build system, to compile the WhiteboxTools project. The `--release` flag optimizes the build for performance, resulting in a smaller and faster executable suitable for distribution or production use. The compiled binary will be located in `target/release/`.

```Shell
cargo build --release
```

--------------------------------

### Building WhiteboxTools Docker Image

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command builds a Docker image for WhiteboxTools. The `-t whitebox-tools` flag tags the image with the name `whitebox-tools`, making it easy to reference. The `-f docker/whitebox-tools.dockerfile` specifies the Dockerfile to use, located within the `docker` subdirectory of the project. The `.` indicates that the build context is the current directory.

```Shell
docker build -t whitebox-tools -f docker/whitebox-tools.dockerfile .
```

--------------------------------

### Printing WhiteboxTools License Information (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

Shows how to access and display the license information associated with the WhiteboxTools library, ensuring compliance and understanding of usage terms.

```Python
print(wbt.license())
```

--------------------------------

### Initializing WhiteboxTools and Setting Working Directory (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This snippet initializes an instance of the WhiteboxTools class and sets its working directory. The working directory is crucial for tools that read or write files, as it defines the default location for input and output data.

```Python
wbt = WhiteboxTools()
wbt.work_dir = "/path/to/data/" # Sets the Whitebox working directory
```

--------------------------------

### Running WhiteboxTools in Docker Container with Data Volume

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command runs the `whitebox-tools` Docker image. The `--rm` flag automatically removes the container after it exits, and `-it` provides an interactive pseudo-TTY. The `-v "/path/to/data/directory/":/data` binds a local data directory to `/data` inside the container, allowing WhiteboxTools to access input files (e.g., `dem.tif`) and write output files (e.g., `out.tif`). The subsequent arguments (`--run=IntegralImage -i=dem.tif -o=out.tif`) are passed directly to the WhiteboxTools executable within the container.

```Shell
docker run --rm -it -v "/path/to/data/directory/":/data whitebox-tools --run=IntegralImage -i=dem.tif -o=out.tif
```

--------------------------------

### Launching WhiteboxTools Runner UI (Python 3)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This command executes the WhiteboxTools Runner, a Tkinter-based user interface, using the `python3` interpreter. It provides a graphical way to interact with WhiteboxTools without direct command-line input.

```Bash
python3 wb_runner.py
```

--------------------------------

### Listing Filtered WhiteboxTools Tools (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

Shows how to list WhiteboxTools tools by applying a filter, in this case, searching for tools with 'lidar' or 'LAS' in their name or description, enabling targeted tool discovery.

```Python
print(wbt.list_tools(['lidar', 'LAS']))
```

--------------------------------

### Calling WhiteboxTools Slope Tool in Python

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This Python snippet illustrates the simplified method for invoking WhiteboxTools functions. It shows the initialization of the `WhiteboxTools` object and then calls the `slope` tool, passing input (`DEM.dep`) and output (`slope.dep`) file paths as arguments. This approach was introduced in version 0.4.0 for improved Python scripting ergonomics.

```Python
wt = WhiteboxTools()
wt.slope(‘DEM.dep’, ‘slope.dep’)
```

--------------------------------

### Changing Directory for Docker Image Build

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command navigates the terminal into the WhiteboxTools project directory. This is a prerequisite for running the `docker build` command, ensuring that the Docker daemon can locate the `docker/whitebox-tools.dockerfile` and the necessary build context.

```Shell
cd /path/to/folder/whitebox-tools/
```

--------------------------------

### Cloning WhiteboxTools Repository for Docker Image Build

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command clones the WhiteboxTools source code repository from GitHub, which is the first step required to build a custom Docker image for WhiteboxTools. The source code is needed by the Docker build process to create the container image.

```Shell
git clone https://github.com/jblindsay/whitebox-tools.git
```

--------------------------------

### Printing WhiteboxTools Version Information (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

Demonstrates how to obtain and print the current version number of the WhiteboxTools library, useful for debugging and ensuring compatibility.

```Python
print("Version information: {}".format(wbt.version()))
```

--------------------------------

### Compiling WhiteboxTools from Source (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/README.md

This command executes the `build.py` Python script, which compiles the WhiteboxTools source code. This script requires a Python environment and may take several minutes to complete, creating a new `WBT` folder containing the compiled program files.

```Python
python build.py
```

--------------------------------

### Launching WhiteboxTools Runner UI (Default Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This command executes the WhiteboxTools Runner, a Tkinter-based user interface, using the default `python` interpreter. It is an alternative to `python3 wb_runner.py` when Python 3 is the system's default Python version.

```Bash
python wb_runner.py
```

--------------------------------

### Executing a WhiteboxTools Function (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This command demonstrates how to call a specific WhiteboxTools function, d_inf_flow_accumulation, from a Python script. It takes input and output file paths as arguments and enables logging for the operation.

```Python
wbt.d_inf_flow_accumulation("DEM.dep", "output.dep", log=True)
```

--------------------------------

### Changing Directory to WhiteboxTools Folder (Shell)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/README.md

This command changes the current working directory in a terminal or command prompt to the `whitebox-tools` folder, which is a prerequisite for compiling the source code. Users must replace `/path/to/folder/` with the actual location where they decompressed the WhiteboxTools source.

```Shell
cd /path/to/folder/whitebox-tools/
```

--------------------------------

### Changing Directory to WhiteboxTools Folder

Source: https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/README.md

This command changes the current working directory in the terminal to the WhiteboxTools project folder. This step is crucial before executing build commands like `cargo build` or `docker build`, as these commands expect to be run from the project's root directory.

```Shell
cd /path/to/folder/whitebox-tools/
```

--------------------------------

### Importing WhiteboxTools Library (Python)

Source: https://github.com/jblindsay/whitebox-tools/blob/master/readme.txt

This line imports the WhiteboxTools class from the whitebox_tools module, making its functionalities available for use within a Python script. It is the first step to programmatically interact with the library.

```Python
from whitebox_tools import WhiteboxTools
```

=== COMPLETE CONTENT === This response contains all available snippets from this library. No additional content exists. Do not make further requests.