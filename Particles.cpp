/* 
 * File:   Particles.cpp
 * Author: gaurav
 * 
 * Created on 4 April, 2013, 9:01 AM
 */

#include "Particles.h"

Particles::Particles(Vec2d& position) {
    
    this->pos[0] = position[0];
    this->pos[1] = position[1];
    
    this->timeAlive = 0.0;
    //this->temp = 
    
}

//Particles::Particles(const Particles& orig) {
//}

Particles::~Particles() {
}

