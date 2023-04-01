# Parallel computing task #2
### Installation of AMD Framewave library
- To use the shared libraries, create the following symbolic links.
For 64 bit installation
```
cd ExampleDir/Framewave/lib
```
For 32 bit installation
```
cd ExampleDir/Framewave/lib
```
Then create the following soft links using the following commands. 
```
ln -sf ./libfwBase.so.1.3.1 libfwBase.so.1
ln -sf ./libfwImage.so.1.3.1 libfwImage.so.1
ln -sf ./libfwJPEG.so.1.3.1 libfwJPEG.so.1
ln -sf ./libfwSignal.so.1.3.1 libfwSignal.so.1
ln -sf ./libfwVideo.so.1.3.1 libfwVideo.so.1
```
```
ln -sf ./libfwBase.so.1 libfwBase.so
ln -sf ./libfwImage.so.1 libfwImage.so
ln -sf ./libfwJPEG.so.1 libfwJPEG.so
ln -sf ./libfwSignal.so.1 libfwSignal.so
ln -sf ./libfwVideo.so.1 libfwVideo.so
```
- Compile a cpp file that uses FW as follows.
```
gcc -O3 -Wall -Werror lab2.c -o lab2-framewave -lfwBase -lfwSignal -LFramewave/lib -lm
```
- Before running the application, make sure the Framewave library is in the environment's shared library (LD_LIBRARY_PATH) search path.
```
export LD_LIBRARY_PATH=ExampleDir/Framewave/lib:$LD_LIBRARY_PATH
```
(example)
```
export LD_LIBRARY_PATH=Framewave/lib:$LD_LIBRARY_PATH
```
