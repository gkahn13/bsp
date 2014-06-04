#ifndef _UTILS_H__
#define _UTILS_H__

#include <unistd.h>
#include <termios.h>

#include <armadillo>
using namespace arma;

namespace utils {

inline char getch() {
	char buf = 0;
	struct termios old = {0};
	if (tcgetattr(0, &old) < 0) {
		perror("tcsetattr()");
	}
	old.c_lflag &= ~ICANON;
	old.c_lflag &= ~ECHO;
	old.c_cc[VMIN] = 1;
	old.c_cc[VTIME] = 0;
	if (tcsetattr(0, TCSANOW, &old) < 0) {
		perror("tcsetattr ICANON");
	}
	if (read(0, &buf, 1) < 0) {
		perror ("read()");
	}
	old.c_lflag |= ICANON;
	old.c_lflag |= ECHO;
	if (tcsetattr(0, TCSADRAIN, &old) < 0) {
		perror ("tcsetattr ~ICANON");
	}
	return (buf);
}

// r,g,b values are from 0 to 1
// h = [0,1], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)

inline mat rgb_to_hsv(const mat &rgb) {
	mat hsv(rgb.n_rows, rgb.n_cols, fill::zeros);
	double min, max, delta;

	min = rgb.min();
	max = rgb.max();
	hsv(2) = max; // v

	delta = max - min;

	if( max != 0 ) {
		hsv(1) = delta / max;		// s
	}
	else {
		// r = g = b = 0		// s = 0, v is undefined
		hsv(1) = 0;
		hsv(0) = -1;

		return hsv;
	}

	if( rgb(0) == max ) {
		hsv(0) = ( rgb(1) - rgb(2) ) / delta;		// between yellow & magenta
	}
	else if( rgb(1) == max ) {
		hsv(0) = 2 + ( rgb(2) - rgb(0) ) / delta;	// between cyan & yellow
	}
	else {
		hsv(0) = 4 + ( rgb(0) - rgb(1) ) / delta;	// between magenta & cyan
	}

	hsv(0) *= 60;				// degrees
	if( hsv(0) < 0 ) {
		hsv(0) += 360;
	}

	hsv(0) /= 360;
	return hsv;
}

inline mat hsv_to_rgb(const mat &hsv) {
	mat rgb(hsv.n_rows, hsv.n_cols, fill::zeros);
	double h = hsv(0), s = hsv(1), v = hsv(2);
	int i;
	double f, p, q, t;

	if( s == 0 ) {
		// achromatic (grey)
		rgb(0) = rgb(1) = rgb(2) = v;
		return rgb;
	}

	h *= 360; // [0,1] --> [0,360]
	h /= 60;  // sector 0 to 5
	i = floor( h );
	f = h - i;  // factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );

	switch( i ) {
		case 0:
			rgb(0) = v;
			rgb(1) = t;
			rgb(2) = p;
			break;
		case 1:
			rgb(0) = q;
			rgb(1) = v;
			rgb(2) = p;
			break;
		case 2:
			rgb(0) = p;
			rgb(1) = v;
			rgb(2) = t;
			break;
		case 3:
			rgb(0) = p;
			rgb(1) = q;
			rgb(2) = v;
			break;
		case 4:
			rgb(0) = t;
			rgb(1) = p;
			rgb(2) = v;
			break;
		default:		// case 5:
			rgb(0) = v;
			rgb(1) = p;
			rgb(2) = q;
			break;
	}

	return rgb;
}

}

#endif
