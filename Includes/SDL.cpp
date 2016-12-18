#include <fstream>
#include "SDL.h"

using namespace std;

SDL::SDL(string path) {
    ifstream ifs;
    ifs.open (path, ifstream::in);
    this->supersampling = false;
    while(!ifs.eof()){
        if(ifs.peek() == '\n' || ifs.peek() == ' ') {
            ifs.get();
        }
        else if(ifs.peek() == '#'){
            ifs.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        else {
            string tag;
            ifs >> tag;
            if(!tag.compare("output")){
                ifs >> this->output;
            }
            else if(!tag.compare(("eye"))){
                T3 eye = T3();
                ifs >> eye.x
                    >> eye.y
                    >> eye.z;
                this->eye = eye;
            }
            else if(!tag.compare("ortho")){
                Ortho ortho = Ortho();
                ifs >> ortho.x0
                    >> ortho.y0
                    >> ortho.x1
                    >> ortho.y1;
                this->ortho = ortho;
            }
            else if(!tag.compare(("size"))){
                Size size = Size();
                ifs >> size.w
                    >> size.h;
                this->size = size;
            }
            else if(!tag.compare("background")){
                T3 background = T3();
                ifs >> background.x
                    >> background.y
                    >> background.z;
                this->background = background;
            }
            else if(!tag.compare("ambient")){
                ifs >> this->ambient;
            }
            else if(!tag.compare("light")){
                Light light = Light();
                ifs >> light.coords.x
                    >> light.coords.y
                    >> light.coords.z;
                this->lights.push_back(light);
            }
            else if(!tag.compare("supersample")) {
                string decisor;
                ifs >> decisor;
                decisor.compare("on") ? this->supersampling = false : this->supersampling = true;
            }
            else if(!tag.compare("profundidade")){
                ifs >> this->depth;
            }
            else if(!tag.compare("object")){
               // Obj o = Obj();
               /* ifs >> o.name
                    >> o.color.r
                    >> o.color.g
                    >> o.color.b
                    >> o.ka
                    >> o.kd
                    >> o.ks
                    >> o.n
                    >> o.KS
                    >> o.KT
                    >> o.ir;*/
                Object object = Object();
                string type;
                ifs >> object.a
                    >> object.b
                    >> object.c
                    >> object.d
                    >> object.e
                    >> object.f
                    >> object.g
                    >> object.h
                    >> object.j
                    >> object.k
                    >> object.color.r
                    >> object.color.g
                    >> object.color.b
                    >> object.ka
                    >> object.kd
                    >> object.ks
                    >> object.n
                    >> object.KS
                    >> object.KT
                    >> object.ir;
                this->objects.push_back(object);
            }
        }
    }
    ifs.close();
}

string SDL::getOutput() {
    return this->output;
}

T3 SDL::getEye(){
    return this->eye;
}

Ortho SDL::getOrtho(){
    return this->ortho;
}

Size SDL::getSize(){
    return this->size;
}

T3 SDL::getBackground(){
    return this->background;
}

double SDL::getAmbient(){
    return this->ambient;
}

vector<Light> SDL::getLights(){
    return this->lights;
}

bool SDL::getSuperSampling(){
    return this->supersampling;
}

double SDL::getDepth(){
    return this->depth;
}

vector<Object> SDL::getObjects(){
    return this->objects;
}
