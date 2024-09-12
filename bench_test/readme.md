# bench_test
## Data path
Download the ASIC data that you need from the google drive and store it on your machine. Let's say the path you choose is `<path-to-your-datafolder>`. Now you have do the following to access the data with your julia code: 
1. Adding environmental variable to `~/.bashrc`:
```
export ASIC_DATA="<path-to-your-datafolder>"  
```
2. Adding environmental variable to VS code settings. Go to `terminal.integrated.env.osx` (or `terminal.integrated.env.linux` depending on your operating system) and add: 
```
  "ASIC_DATA":"<path-to-your-datafolder>" 
```