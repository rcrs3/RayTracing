#include <bits/stdc++.h>
#include "Includes/header.h"

using namespace std;


int main() {
	SDL* sdl = new SDL("onesphere.sdl");
	cout << sdl->getOutput() << endl;

	delete sdl;
	return 0;
}