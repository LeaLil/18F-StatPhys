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

    output.printInitialTemperature(properties[1] / fac);
    output.printIterationStart();

    performStep(x, v, 0, par.initialTime); //Init stuff...
    /* dynamics step */
    double time = par.initialTime;
    for (int nstep = 0; nstep < par.numberMDSteps; nstep++) {
        time += par.timeStep;
        //performStep(x, v, nstep, time);
        performMetropolisalgorithm(x, v, nstep, time);
    }

    printAverages(time);
}

void MDRun::performMetropolisalgorithm(std::vector<double> &positions, std::vector<double> &velocities, int nstep,

                                       double time) {

    //Select random atom to displace.
    // ExTODO: Is that really what we're supposed to do? => Yes, according to "Computer simulations for liquids" by Michael P. Allen
    int nrOfAtoms = positions.size() / 3;
    int atomAtXPosition = (rand() % nrOfAtoms) * 3;
    double coords[3] = {positions.at(atomAtXPosition), positions.at(atomAtXPosition + 1),
                        positions.at(atomAtXPosition + 2)};
    double currentKineticEnergy = properties[1];
    //Compute current potential energy
    forceCalculator.calculate(positions, forces);
    double potentialEnergyOldSystem = forceCalculator.getPotentialEnergy();

    //For metropolis we have to compute y = x_i + r*q. R is supposedly the size of the box
    bool noAcceptedNewConfiguration = true;
    //TODO: Choose r s.t. acceptance rate about 0.5

    while (noAcceptedNewConfiguration) {
        //random Number between -1 and 1
        for (int i = 0; i < 3; i++) {
            positions.at(atomAtXPosition + i) = positions.at(atomAtXPosition + i) + r * getRandomNumberForMetropolis();
        }

        //Compute potentialEnergy of new system
        PeriodicBoundaryConditions::recenterAtoms(par.numberAtoms, positions, par.boxSize);

        forceCalculator.calculate(positions, forces);
        double potentialEnergyNewSystem = forceCalculator.getPotentialEnergy();
        double deltaPotentialEnergy = potentialEnergyNewSystem - potentialEnergyOldSystem;

        //exTODO: What's the behavior if T = 0? -> Obviously nothing is moving at T = 0
        double probabilityOfGettingAccepted = std::exp(
                -deltaPotentialEnergy / (boltzmannConstant * (par.targetTemperature))); //WE assume that temperature is constant, so shouldn't be an issue.

        bool accepted = probabilityOfGettingAccepted > ((double) rand() / (RAND_MAX)) || deltaPotentialEnergy < 0;

        if (accepted) {
            nrOfAcceptedConfigurations++;
            noAcceptedNewConfiguration = false;
            properties[2] = properties[2] - potentialEnergyOldSystem + potentialEnergyNewSystem; //potential Energy
            radialDistribution.addInstantaneousDistribution(forceCalculator.getInstantaneousRadialDistribution());
            properties[3] = forceCalculator.getVirial();
            properties[5] = computeVelocityScalingFactor();

            double oldKineticEnergy;
            double newKineticEnergy;
            performLeapFrog(properties[5], positions, velocities, oldKineticEnergy, newKineticEnergy);

            properties[1] = newKineticEnergy; //New kinetic energy is total energy - new potential energy
            properties[4] = 2. * (properties[1] - properties[3]) / (vol * 3.);
            properties[0] = properties[1] + properties[2];


            if (par.mdType == SimulationType::constantTemperature) {
                ekg = currentKineticEnergy;
            }

        } else {
            nrOfRejectedConfigurations++;
            //Reset positions for new run
            for (int i = 0; i < 3; i++) {
                positions.at(atomAtXPosition + i) = coords[i];
            }

        }
    }
    double ratioOfAcceptedVsNotAccepted =
            nrOfAcceptedConfigurations / (nrOfAcceptedConfigurations + nrOfRejectedConfigurations);
    if (ratioOfAcceptedVsNotAccepted > 0.55) {
        r *= 1.05;
    }else if (ratioOfAcceptedVsNotAccepted < 0.45) {
        r *= 0.95;
    }

    /* update arrays for averages and fluctuations */
    for (int m = 0; m < numberProperties; m++) {
        averages[m] += properties[m];
        fluctuations[m] += properties[m] * properties[m];
    }

    printOutputForStep(positions, velocities, nstep, time);

    //Compute stuff in new configuration that we haven't computed yet
    //performStep(positions, velocities,nstep, time);

}

float MDRun::getRandomNumberForMetropolis() const { return -1 + static_cast <double> (rand()) /
                                                                (static_cast <double> (RAND_MAX / (2)));
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

    printOutputForStep(positions, velocities, nstep, time);
}

void MDRun::performLeapFrog(double scal, std::vector<double> &positions, std::vector<double> &velocities,
                            double &oldKineticEnergy, double &newKineticEnergy) const {
    oldKineticEnergy= 0.;
    newKineticEnergy= 0.;/* perform leap-frog integration step,
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
