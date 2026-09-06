//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# This program is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//#############################################################################

#include "remap_plan_internal.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace Loci {
  namespace detail {

    bool buildDistributedCellRemapPlan(CPTR<AMRRemapPlan>& plan,
          AMRRemapReport& report, gatherCommSchedule& sourceGather,
          const std::vector<std::pair<int, int>>& targetSource,
          const store<double>& sourceVolume,
          const_store<vector3d<double>>& sourceCenter,
          const_store<CellId>& sourceCellId,
          const store<int>& sourceTargetCount, dataPartitionP sourcePartition,
          const entitySet& localTargetCells, const Map& targetLocalToGlobal,
          const const_store<vector3d<double>>& targetCenter,
          const const_store<double>& targetVolume,
          const const_store<CellId>& targetCellId,
          const store<int>& refinedSource,
          const multiStore<int>& refinedSourceToTarget,
          const store<double>& gatheredTargetVolume,
          const store<vector3d<double>>& gatheredTargetCenter, MPI_Comm comm) {
      const double tolerance = 1.0e-8 ;
      plan = static_cast<AMRRemapPlan*>(0) ;
      report = AMRRemapReport() ;

      int localIdentityDomainsValid =
            sourceCellId.domain() == sourceVolume.domain() &&
                        (localTargetCells - targetCellId.domain()).size() == 0
                  ? 1
                  : 0 ;
      int identityDomainsValid = 0 ;
      MPI_Allreduce(&localIdentityDomainsValid, &identityDomainsValid, 1,
            MPI_INT, MPI_MIN, comm) ;
      if (identityDomainsValid == 0) {
        report.invalidIdentities = 1 ;
        return false ;
      }

      entitySet requestedSources ;
      for(size_t i=0;i<targetSource.size();++i)
        requestedSources += targetSource[i].second ;
      sourceGather.generateSchedule(requestedSources,sourcePartition) ;

      store<vector3d<double> > reconstructionPoint ;
      reconstructionPoint.allocate(sourceVolume.domain()) ;
      FORALL(sourceVolume.domain(),source) {
        reconstructionPoint[source] = sourceCenter[source] ;
      } ENDFORALL ;
      FORALL(refinedSource.domain(),source) {
        double totalVolume = 0.0 ;
        vector3d<double> weightedCenter(0.0,0.0,0.0) ;
        for(int target=0;
            target<refinedSourceToTarget.vec_size(source);++target) {
          const int localTarget = refinedSourceToTarget[source][target] ;
          totalVolume += gatheredTargetVolume[localTarget] ;
          weightedCenter += gatheredTargetVolume[localTarget]*
            gatheredTargetCenter[localTarget] ;
        }
        if(totalVolume > 0.0)
          reconstructionPoint[refinedSource[source]] =
            weightedCenter/totalVolume ;
      } ENDFORALL ;

      store<double> importedSourceVolume ;
      store<vector3d<double> > importedSourceCenter ;
      store<vector3d<double> > importedReconstructionPoint ;
      store<CellId> importedSourceCellId ;
      store<int> importedSourceTargetCount ;
      sourceGather.gatherData(importedSourceVolume,sourceVolume) ;
      sourceGather.gatherData(importedSourceCenter,sourceCenter) ;
      sourceGather.gatherData(importedReconstructionPoint,
                              reconstructionPoint) ;
      sourceGather.gatherData(importedSourceCellId, sourceCellId) ;
      sourceGather.gatherData(importedSourceTargetCount,
                              sourceTargetCount) ;

      std::map<int,int> sourceToLocal ;
      Loci::getLocalContextMap(sourceToLocal,requestedSources) ;
      std::map<int,Entity> targetToLocal ;
      FORALL(localTargetCells,target) {
        targetToLocal[targetLocalToGlobal[target]] = target ;
      } ENDFORALL ;

      std::map<int,size_t> targetDegree ;
      for(size_t i=0;i<targetSource.size();++i)
        targetDegree[targetSource[i].first]++ ;

      std::vector<AMRCellGeometry> sourceGeometry ;
      std::vector<size_t> sourceDegrees ;
      sourceGeometry.reserve(requestedSources.size()) ;
      sourceDegrees.reserve(requestedSources.size()) ;
      for(entitySet::const_iterator source=requestedSources.begin();
          source!=requestedSources.end();++source) {
        const int localSource = sourceToLocal[*source] ;
        sourceGeometry.push_back(
              AMRCellGeometry(*source, importedSourceCellId[localSource],
                    importedSourceVolume[localSource],
                    importedSourceCenter[localSource],
                    importedReconstructionPoint[localSource])) ;
        sourceDegrees.push_back(
          size_t(importedSourceTargetCount[localSource])) ;
      }

      std::vector<AMRCellGeometry> targetGeometry ;
      targetGeometry.reserve(targetToLocal.size()) ;
      for(std::map<int,Entity>::const_iterator target=targetToLocal.begin();
          target!=targetToLocal.end();++target)
        targetGeometry.push_back(AMRCellGeometry(target->first,
              targetCellId[target->second], targetVolume[target->second],
              targetCenter[target->second])) ;

      std::vector<AMRCellContribution> contributions ;
      contributions.reserve(targetSource.size()) ;
      for(size_t i=0;i<targetSource.size();++i) {
        const int target = targetSource[i].first ;
        const int source = targetSource[i].second ;
        const int localSource = sourceToLocal[source] ;
        const std::map<int,Entity>::const_iterator localTarget =
          targetToLocal.find(target) ;
        if(localTarget == targetToLocal.end())
          continue ;

        double overlapVolume = targetVolume[localTarget->second] ;
        vector3d<double> overlapCentroid = targetCenter[localTarget->second] ;
        if(targetDegree[target] > 1) {
          overlapVolume = importedSourceVolume[localSource] ;
          overlapCentroid = importedSourceCenter[localSource] ;
        }
        contributions.push_back(AMRCellContribution(source, target,
              importedSourceCellId[localSource],
              targetCellId[localTarget->second], overlapVolume,
              overlapCentroid)) ;
      }

      AMRRemapReport localReport ;
      CPTR<AMRRemapPlan> localPlan = AMRRemapPlan::createCellPlanPartition(
        sourceGeometry,sourceDegrees,targetGeometry,contributions,
        localReport,tolerance) ;

      size_t sourceCoverageErrors = 0 ;
      size_t sourceMomentErrors = 0 ;
      double maximumSourceVolumeError = 0.0 ;
      double maximumSourceCentroidError = 0.0 ;
      FORALL(sourceTargetCount.domain(),source) {
        if(sourceTargetCount[source] == 0)
          sourceCoverageErrors++ ;
      } ENDFORALL ;

      // A one-target source occurs on exactly one target owner. Refined
      // sources are checked once on their source owner from the complete
      // gathered target set.
      for(size_t source=0;source<sourceGeometry.size();++source) {
        if(sourceDegrees[source] != 1)
          continue ;
        for(size_t entry=0;entry<contributions.size();++entry) {
          if(contributions[entry].sourceCell != sourceGeometry[source].cell)
            continue ;
          const double volumeError = std::abs(
            contributions[entry].overlapVolume-sourceGeometry[source].volume) ;
          if(volumeError > tolerance*
             std::max(1.0,std::abs(sourceGeometry[source].volume)))
            sourceCoverageErrors++ ;
          maximumSourceVolumeError =
            std::max(maximumSourceVolumeError,volumeError) ;

          const double centroidError = norm(
            contributions[entry].overlapCentroid-
            sourceGeometry[source].reconstructionPoint) ;
          if(centroidError > tolerance*std::max(
               1.0,norm(sourceGeometry[source].reconstructionPoint)))
            sourceMomentErrors++ ;
          maximumSourceCentroidError =
            std::max(maximumSourceCentroidError,centroidError) ;
          break ;
        }
      }
      FORALL(refinedSource.domain(),source) {
        double coveredVolume = 0.0 ;
        for(int target=0;
            target<refinedSourceToTarget.vec_size(source);++target)
          coveredVolume += gatheredTargetVolume[
            refinedSourceToTarget[source][target]] ;
        const double volumeError =
          std::abs(coveredVolume-sourceVolume[refinedSource[source]]) ;
        if(volumeError > tolerance*std::max(
             1.0,std::abs(sourceVolume[refinedSource[source]])))
          sourceCoverageErrors++ ;
        maximumSourceVolumeError =
          std::max(maximumSourceVolumeError,volumeError) ;
      } ENDFORALL ;

      unsigned long long localCounts[11] = {
            static_cast<unsigned long long>(sourceVolume.domain().size()),
            static_cast<unsigned long long>(localTargetCells.size()),
            static_cast<unsigned long long>(targetSource.size()),
            static_cast<unsigned long long>(localReport.invalidGeometry),
            static_cast<unsigned long long>(localReport.invalidIdentities),
            static_cast<unsigned long long>(localReport.duplicateContributions),
            static_cast<unsigned long long>(sourceCoverageErrors),
            static_cast<unsigned long long>(localReport.missingTargetCells),
            static_cast<unsigned long long>(sourceMomentErrors),
            static_cast<unsigned long long>(
                  localReport.inconsistentTargetMoments),
            static_cast<unsigned long long>(localReport.unsupportedRelations)} ;
      unsigned long long globalCounts[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} ;
      MPI_Allreduce(localCounts, globalCounts, 11, MPI_UNSIGNED_LONG_LONG,
            MPI_SUM, comm) ;

      double localErrors[4] = {
        maximumSourceVolumeError,
        localReport.maximumTargetVolumeError,
        maximumSourceCentroidError,
        localReport.maximumTargetCentroidError
      } ;
      double globalErrors[4] = {0.0,0.0,0.0,0.0} ;
      MPI_Allreduce(localErrors,globalErrors,4,MPI_DOUBLE,MPI_MAX,comm) ;

      int localPlanValid =
        localPlan != static_cast<AMRRemapPlan*>(0) ? 1 : 0 ;
      int allPlansValid = 0 ;
      MPI_Allreduce(&localPlanValid,&allPlansValid,1,MPI_INT,MPI_MIN,comm) ;

      report.sourceCells = size_t(globalCounts[0]) ;
      report.targetCells = size_t(globalCounts[1]) ;
      report.contributions = size_t(globalCounts[2]) ;
      report.invalidGeometry = size_t(globalCounts[3]) ;
      report.invalidIdentities = size_t(globalCounts[4]) ;
      report.duplicateContributions = size_t(globalCounts[5]) ;
      report.missingSourceCells = size_t(globalCounts[6]) ;
      report.missingTargetCells = size_t(globalCounts[7]) ;
      report.inconsistentSourceMoments = size_t(globalCounts[8]) ;
      report.inconsistentTargetMoments = size_t(globalCounts[9]) ;
      report.unsupportedRelations = size_t(globalCounts[10]) ;
      report.sourceCoverageChecked = true ;
      report.maximumSourceVolumeError = globalErrors[0] ;
      report.maximumTargetVolumeError = globalErrors[1] ;
      report.maximumSourceCentroidError = globalErrors[2] ;
      report.maximumTargetCentroidError = globalErrors[3] ;
      report.valid = allPlansValid != 0 && report.invalidGeometry == 0 &&
                     report.invalidIdentities == 0 &&
                     report.duplicateContributions == 0 &&
                     report.missingSourceCells == 0 &&
                     report.missingTargetCells == 0 &&
                     report.inconsistentSourceMoments == 0 &&
                     report.inconsistentTargetMoments == 0 &&
                     report.unsupportedRelations == 0 ;
      if(report.valid)
        plan = localPlan ;
      return report.valid ;
    }
  }
}
