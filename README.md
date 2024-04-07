# Build  
cmake -D CMAKE_BUILD_TYPE=Release -DEDO=ON -B build -DENABLE_CMAKE_EXAMPLE=ON

I dir "masteroppgave" kjør:
     cmake --build build
     ./build/hillclimber

For å få tilgang til html-filene:
     gå til mappen "paradiseo/eo/tutorial/html
     cmd: python3 -m http.server 8000