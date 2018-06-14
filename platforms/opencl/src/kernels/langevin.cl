enum {VelScale, ForceScale, NoiseScale, MaxParams};

/**
 * Perform the first step of Langevin integration.
 */

__kernel void integrateLangevinPart1(__global mixed4* restrict velm, __global const real4* restrict force, __global mixed4* restrict posDelta,
        __global const mixed* restrict paramBuffer, __global const mixed2* restrict dt, __global const float4* restrict random, unsigned int randomIndex) {
    mixed vscale = paramBuffer[VelScale];
    mixed fscale = paramBuffer[ForceScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed stepSize = dt[0].y;
    int index = get_global_id(0);
    randomIndex += index;
    while (index < NUM_ATOMS) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed sqrtInvMass = sqrt(velocity.w);
            velocity.x = vscale*velocity.x + fscale*velocity.w*force[index].x + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + fscale*velocity.w*force[index].y + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + fscale*velocity.w*force[index].z + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
            posDelta[index] = stepSize*velocity;
        }
        randomIndex += get_global_size(0);
        index += get_global_size(0);
    }
}


/**
 * Perform the second step of Langevin integration.
 */

__kernel void integrateLangevinPart2(__global real4* restrict posq, __global real4* restrict posqCorrection, __global const mixed4* restrict posDelta, __global mixed4* restrict velm, __global const mixed2* restrict dt) {
#ifdef SUPPORTS_DOUBLE_PRECISION
    double invStepSize = 1.0/dt[0].y;
#else
    float invStepSize = 1.0f/dt[0].y;
    float correction = (1.0f-invStepSize*dt[0].y)/dt[0].y;
#endif
    int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        mixed4 vel = velm[index];
        if (vel.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            mixed4 delta = posDelta[index];
            pos.xyz += delta.xyz;
#ifdef SUPPORTS_DOUBLE_PRECISION
            vel.xyz = convert_mixed4(invStepSize*convert_double4(delta)).xyz;
#else
            vel.xyz = invStepSize*delta.xyz + correction*delta.xyz;
#endif

#ifdef USE_MIXED_PRECISION
            posq[index] = convert_real4(pos);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = vel;
        }
        index += get_global_size(0);
    }
}

/**
 * Select the step size to use for the next step.
 */

__kernel void selectLangevinStepSize(mixed maxStepSize, mixed errorTol, mixed tau, mixed kT, __global mixed2* restrict dt,
        __global const mixed4* restrict velm, __global const real4* restrict force, __global mixed* restrict paramBuffer, __local mixed* restrict params, __local mixed* restrict error) {
    // Calculate the error.

    mixed err = 0.0f;
    unsigned int index = get_local_id(0);
    while (index < NUM_ATOMS) {
        real4 f = force[index];
        mixed invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass*invMass;
        index += get_global_size(0);
    }
    error[get_local_id(0)] = err;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Sum the errors from all threads.

    for (unsigned int offset = 1; offset < get_local_size(0); offset *= 2) {
        if (get_local_id(0)+offset < get_local_size(0) && (get_local_id(0)&(2*offset-1)) == 0)
            error[get_local_id(0)] += error[get_local_id(0)+offset];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_global_id(0) == 0) {
        // Select the new step size.

        mixed totalError = sqrt(error[0]/(NUM_ATOMS*3));
        mixed newStepSize = sqrt(errorTol/totalError);
        mixed oldStepSize = dt[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        dt[0].y = newStepSize;

        // Recalculate the integration parameters.

        mixed vscale = exp(-newStepSize/tau);
        mixed fscale = (1-vscale)*tau;
        mixed noisescale = sqrt(2*kT/tau)*sqrt(0.5f*(1-vscale*vscale)*tau);
        params[VelScale] = vscale;
        params[ForceScale] = fscale;
        params[NoiseScale] = noisescale;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if (get_local_id(0) < MaxParams)
        paramBuffer[get_local_id(0)] = params[get_local_id(0)];
}




//lambda-hybridize forces
__kernel void bedamForce(__global const real4* restrict posq, __global real4* restrict force, real lambdaId){

    int index = get_global_id(0);
    int halfAtoms = NUM_ATOMS/2;
    while (index < halfAtoms) {

        real4 force1 = force[index];
        real4 force2 = force[index+halfAtoms];

        force1.x = lambdaId*force1.x + (1.0 - lambdaId)*force2.x;
        force1.y = lambdaId*force1.y + (1.0 - lambdaId)*force2.y;
        force1.z = lambdaId*force1.z + (1.0 - lambdaId)*force2.z;

        force[index] = force1;

        index += get_global_size(0);
	
    }

}


/**
 * restraint force between two groups of atoms (Vsite)
 */
__kernel void CMForceCalcKernel(__global const int*  restrict particle_indexes1,
				__global const int*  restrict particle_indexes2,
				__global const real4* restrict posq, //atomic positions
				__global real4* restrict force,
				int numParticles1, int numParticles2,
				real kf,
				real r0){
				
  uint id = get_global_id(0);
  if(id == 0 && numParticles1 > 0 && numParticles2 > 0){
    real4 rf1 = (real4)(0,0,0,0);
    real4 rf2 = (real4)(0,0,0,0);
    
    real4 cm1 = (real4)(0,0,0,0);
    float vn1 = 1./numParticles1;
    for ( int i = 0; i < numParticles1; i++ ){
      int index = particle_indexes1[i];
      cm1 += vn1*posq[index];
    }
    real4 cm2 = (real4)(0,0,0,0);
    float vn2 = 1./numParticles2;
    for ( int i = 0; i < numParticles2; i++ ){
      int index = particle_indexes2[i];
      cm2 += vn2*posq[index];
    }
    real4 dist = cm2 - cm1;
    float dr = sqrt(dot(dist,dist));
    if(dr>r0) 
    {
      float a = kf*(dr-r0)/dr;
      rf1 = dist *  vn1*a;
      rf2 = dist * -vn2*a;
    }
    
    for ( int i = 0; i < numParticles1; i++ ){
      int index = particle_indexes1[i];
      force[index] += rf1;
    }
    for ( int i = 0; i < numParticles2; i++ ){
      int index = particle_indexes2[i];
      force[index] += rf2;
    }
    
    
  }
}



/**
 * copy the original positions and velocities to the imaginary parts.
 */
__kernel void copyDataToSecondPart(__global real4* restrict posq, __global real4* restrict velm, int ligId) {

    
    //int index = NUM_ATOMS/2;
    int index = get_global_id(0);
    int halfAtoms = NUM_ATOMS/2;
    int ligN =  halfAtoms+ligId;
    
    while (index < NUM_ATOMS) {
	 
	if (index >=halfAtoms) {

	real4 pos1 = posq[index];
	real4 pos2 = posq[index - halfAtoms];
	real4 vel1 = velm[index];
	real4 vel2 = velm[index - halfAtoms];
	pos1.x = (index < ligN)? (pos2.x + 200.0f) : (pos2.x + 100.0f);
 	pos1.y = pos2.y;
	pos1.z = pos2.z;
	vel1.xyz = vel2.xyz;
	posq[index] = pos1;
	velm[index] = vel1;
        }
        //index++;
	index += get_global_size(0);
    }
   



}				   

