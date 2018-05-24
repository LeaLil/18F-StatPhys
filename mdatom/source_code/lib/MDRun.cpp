#include "MDRun.h"
#include "PeriodicBoundaryConditions.h"
#include "CenterOfMassCalculator.h"
#include "TrajectoryFileWriter.h"
#include <cmath>
#include <random>
#include <algorithm>

MDRun::MDRun(const MDParameters &parameters, MDRunOutput &out, TrajectoryFileWriter &trajectoryFileWriter)
        : par(parameters),
          output(out),
          trajectoryWriter(trajectoryFileWriter),
          forceCalculator(parameters),
          radialDistribution(parameters.numberRadialDistrPoints, parameters.radialDistrCutoffRadius) {
}

void MDRun::run(std::vector<double> &x, std::vector<double> &v) {
    forces.resize(x.size());
    synchronizedPositions.resize(x.size());
    radialDistribution.setZero();

    initializeVariables();
    initializeTemperature(v);
    //output.printInitialTemperature(properties[1] / fac);
    //output.printIterationStart();

    bool MConeParticleOnly = true; // This is a hack; needs proper parametrisation

    //Recenter atoms if necessary (e.g. if coords.inp has coordinates on the outside of our box)
    //Not necessary for Leap-Frog, but whatever :)
    PeriodicBoundaryConditions::recenterAtoms(par.numberAtoms, x, par.boxSize);
    forceCalculator.calculate(x, v);
    properties[2] = forceCalculator.getPotentialEnergy();



    /* dynamics step */
    double time = par.initialTime;
    //std::cout << "Metropolis" << std::endl;
    //std::cout << "Normal" << std::endl;
    for (int nstep = 0; nstep < par.numberMDSteps; nstep++) {
        time += par.timeStep;
        if (par.isMonteCarlo) {
            performMetropolisStep(x, v, nstep, time, MConeParticleOnly);
        } else {
            performStep(x, v, nstep, time);

        }

    }

    //  std::cout << "Accepted: " << nrOfAcceptedConfigurations << " Not: " << nrOfRejectedConfigurations;

//printAverages(time);
}

void
MDRun::performMetropolisStep(std::vector<double> &positions, std::vector<double> &velocities, int nstep, double time,
                             bool MConeParticleOnly) {


    int atomPositionStart = 0;
    int atomPositionRange = positions.size();

    // Restrict range for changing positions to 1 atom if (oneParticleOnly)
    // otherwise we'll change positions of all (presumably N) atoms
    if (MConeParticleOnly) {
        //Select random atom to displace.
        // ExTODO: Is that really what we're supposed to do? => Yes, according to "Computer simulations for liquids" by Michael P. Allen
        atomPositionStart = (rand() % par.numberAtoms) * 3;
        atomPositionRange = atomPositionStart + 3;
        // double coords[3] = {positions.at(atomAtXPosition),
        //                     positions.at(atomAtXPosition + 1),
        //                     positions.at(atomAtXPosition + 2)};
    }

    //current potential energy
    double potentialEnergyPreMCStep = properties[2];

    //Make a copy of pre MC step positions in case MC step gets rejected
    std::vector<double> positionsPreMCStep(positions);


    //Get ready for Metropolis MC step
    bool acceptNewConfiguration = false;
    double deltaPotentialEnergy;

    //Random number generator with two rvars
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> testProb(0, 1);
    std::uniform_real_distribution<double> stepScale(-1, 1);

    //For metropolis we have to compute y = x_i + delta_r * step.
    //random Number between -1 and 1
    for (int i = atomPositionStart; i < atomPositionRange; i++) {
        positions.at(i) = positions.at(i) + par.delta_r * stepScale(gen);
    }

    //Compute potential energy post MC step
    PeriodicBoundaryConditions::recenterAtoms(par.numberAtoms, positions, par.boxSize);
    forceCalculator.calculate(positions, forces);
    double potentialEnergyPostMCStep = forceCalculator.getPotentialEnergy();

    //Compute energy difference
    deltaPotentialEnergy = potentialEnergyPostMCStep - potentialEnergyPreMCStep;
    //If energy has been reduced: Smooth sailing accepting step
    if (deltaPotentialEnergy <= 0) {
        acceptNewConfiguration = true;
    }
        //Hmm, energy has been increased: Play the casino
    else {
        //MC Metropolis test
        //exTODO: What's the behavior if T = 0? -> Obviously nothing is moving at T = 0
        double MCtest = std::exp(
                -deltaPotentialEnergy / (boltzmannConstant *
                                         par.targetTemperature)); //WE assume that temperature is constant, so shouldn't be an issue.
        acceptNewConfiguration = (MCtest > testProb(gen));
    }
    //We have an accepted configuration: Increase count and we're done
    if (acceptNewConfiguration) {
        nrOfAcceptedConfigurations++;
        properties[2] = potentialEnergyPostMCStep;
    }
        //This step didn't yield an accepted configuration: Increase count and restore positions
    else {
        nrOfRejectedConfigurations++; // not needed if run as step in context of MDRun::run
        positions = positionsPreMCStep;

    }

    //Calculate running acceptance ratio so delta_r can be tuned for future runs
    double ratioOfAcceptedVsNotAccepted =
            nrOfAcceptedConfigurations / (nrOfAcceptedConfigurations + nrOfRejectedConfigurations);

    // if (ratioOfAcceptedVsNotAccepted > 0.51) {
    //     par.delta_r *= 1.01;
    // } else if (ratioOfAcceptedVsNotAccepted < 0.49) {
    //     par.delta_r *= 0.99;
    // }

    // std::cout << "delta_r for next step : " << par.delta_r << std::endl;


    if (!par.showDistributionInsteadOfCSV) {
        std::cout << nstep << "," << properties[2] << "," << deltaPotentialEnergy << "," << par.delta_r << ","
                  << ratioOfAcceptedVsNotAccepted << std::endl;

    }
}


void MDRun::initializeVariables() {
    nat3 = 3 * par.numberAtoms;

    fac = nat3 * boltzmannConstant / 2.;
    ekin0 = fac * par.targetTemperature;
    halfTimeStep = par.timeStep / 2;
    dtm = par.timeStep / par.atomicMass;
    vol = par.boxSize[0] * par.boxSize[1] * par.boxSize[2];

    nhpr = 100 * par.propertyPrintingInterval;
    nlsq = par.numberMDSteps / 10;
    if (nlsq < 10) {
        nlsq = 10;
    }
    if (nlsq < par.propertyPrintingInterval) {
        nlsq = par.propertyPrintingInterval;
    }
    if (nhpr > nlsq) {
        nhpr = nlsq;
    }
    for (int i = 0; i < numberProperties; i++) {
        properties[i] = 0.;
        averages[i] = 0.;
        fluctuations[i] = 0.;
    }
}

void MDRun::initializeTemperature(const std::vector<double> &velocities) {
    double kineticEnergy = 0;
    for (int j3 = 0; j3 < nat3; j3++) {
        kineticEnergy += velocities[j3] * velocities[j3];
    }
    kineticEnergy *= (par.atomicMass / 2.);
    properties[1] = kineticEnergy;
    if (par.mdType == SimulationType::constantTemperature) {
        if (kineticEnergy < 1.e-6) {
            ekg = ekin0;
        } else {
            ekg = kineticEnergy;
        }
    }
}

void MDRun::performStep(std::vector<double> &positions, std::vector<double> &velocities, int nstep, double time) {
    /* put atoms in central periodic box */

    PeriodicBoundaryConditions::recenterAtoms(par.numberAtoms, positions, par.boxSize);
    double currentPotEnergy = forceCalculator.getPotentialEnergy();
    /* calculate forces, potential energy, virial
     * and contribution to the radial distribution function
     */
    forceCalculator.calculate(positions, forces);
    radialDistribution.addInstantaneousDistribution(forceCalculator.getInstantaneousRadialDistribution());
    double vir = forceCalculator.getVirial();
    properties[2] = forceCalculator.getPotentialEnergy();
    properties[3] = vir;

    /* determine velocity scaling factor, when coupling to a bath */
    double scal = 1;
    if (par.mdType == SimulationType::constantTemperature) {
        double dtt = par.timeStep / par.temperatureCouplingTime;
        scal = std::sqrt(1 + dtt * (ekin0 / ekg - 1));
    }
    double oldKineticEnergy;
    double newKineticEnergy;
    performLeapFrog(scal, positions, velocities, oldKineticEnergy, newKineticEnergy);

    properties[1] = newKineticEnergy;
    properties[0] = properties[1] + properties[2];
    double pres = 2. * (newKineticEnergy - vir) / (vol * 3.);
    properties[4] = pres;
    properties[5] = scal;
    if (par.mdType == SimulationType::constantTemperature) {
        ekg = oldKineticEnergy;
    }

    /* update arrays for averages and fluctuations */
    for (int m = 0; m < numberProperties; m++) {
        averages[m] += properties[m];
        fluctuations[m] += properties[m] * properties[m];
    }

    //printOutputForStep(positions, velocities, nstep, time)
    if (!par.showDistributionInsteadOfCSV) {
        std::cout << nstep << "," << properties[2] << "," << (properties[2] - currentPotEnergy) << std::endl;
    }

}

void MDRun::performLeapFrog(double scal, std::vector<double> &positions, std::vector<double> &velocities,
                            double &oldKineticEnergy, double &newKineticEnergy) const {
    oldKineticEnergy = 0.;
    newKineticEnergy = 0.;/* perform leap-frog integration step,
     * calculate kinetic energy at time t-dt/2 and at time t,
     * and calculate pressure
     */for (int j3 = 0; j3 < nat3; j3++) {
        double oldVelocity = velocities[j3];
        double newVelocity = (oldVelocity + forces[j3] * dtm) * scal;
        oldKineticEnergy += newVelocity * newVelocity;
        newKineticEnergy += (oldVelocity + newVelocity) * (oldVelocity + newVelocity);
        velocities[j3] = newVelocity;
        positions[j3] += newVelocity * par.timeStep;
    }
    oldKineticEnergy *= (par.atomicMass / 2.);
    newKineticEnergy *= (par.atomicMass / 8.);
}

double MDRun::computeVelocityScalingFactor() const {/* determine velocity scaling factor, when coupling to a bath */
    double scal = 1;
    if (par.mdType == SimulationType::constantTemperature) {
        double dtt = par.timeStep / par.temperatureCouplingTime;
        scal = sqrt(1 + dtt * (ekin0 / ekg - 1));
    }
    return scal;
}


void MDRun::printOutputForStep(const std::vector<double> &positions, const std::vector<double> &velocities, int nstep,
                               double time) {
    if ((nstep + 1) == (nstep + 1) / par.trajectoryOutputInterval * par.trajectoryOutputInterval) {
        trajectoryWriter.writeOutTrajectoryStep(positions);
    }

    if (nstep == (nstep + 1) / nhpr * nhpr) {
        output.printPropertiesHeader();
    }

    if ((nstep + 1) == (nstep + 1) / par.propertyPrintingInterval * par.propertyPrintingInterval || nstep == 0) {
        output.printProperties(nstep, time, properties);
    }


    /* calculate and print center of mass motion
     * once in nlsq steps, at time t-dt/2
     * The positions must be back-calculated for t-dt/2, because of the time shift between x and v (leap-frog)
     */
    if ((nstep + 1) == (nstep + 1) / nlsq * nlsq) {
        for (int j3 = 0; j3 < nat3; j3++) {
            synchronizedPositions[j3] = positions[j3] - velocities[j3] * halfTimeStep;
        }
        CenterOfMassCalculator cm;
        cm.update(par.numberAtoms, synchronizedPositions, velocities, par.atomicMass);
        cm.printResults(output);
    }
}

void MDRun::printAverages(double time) {
    double tspan = par.numberMDSteps;
    for (int m = 0; m < numberProperties; m++) {
        averages[m] = averages[m] / tspan;
        fluctuations[m] = std::sqrt(std::abs(fluctuations[m] / tspan - averages[m] * averages[m]));
    }
    output.printAverages(par.numberMDSteps, time, averages);
    output.printRMSFluctuations(par.numberMDSteps, time, fluctuations);
    output.printAverageAndRMSTemperature(averages[1] / fac, fluctuations[1] / fac);
}

const AveragedRadialDistribution &MDRun::getRadialDistribution() const {
    return radialDistribution;
}
