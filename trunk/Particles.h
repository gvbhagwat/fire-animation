/* 
 * File:   Particles.h
 * Author: gaurav
 *
 * Created on 4 April, 2013, 9:01 AM
 */

#ifndef PARTICLES_H
#define	PARTICLES_H

#include "Vec.h"

class Particles {
public:
    Particles(Vec2d&);
    //Particles(const Particles& orig);
    virtual ~Particles();
    
    Vec2d pos;
    double timeAlive;
    double y;
    double temp;


};

#endif	/* PARTICLES_H */

