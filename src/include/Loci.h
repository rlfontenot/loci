#ifndef LOCI_H
#define LOCI_H

#include <Map.h>
#include <store.h>
#include <storeVec.h>
#include <parameter.h>
#include <constraint.h>
#include <path.h>
#include <rule.h>
#include <options_list.h>
#include <fpe.h>
#include <scheduler.h>

using Loci::entitySet ;
using Loci::EMPTY ;
using Loci::interval ;

using Loci::Map ;
using Loci::const_Map ;
using Loci::MapVec ;
using Loci::const_MapVec ;
using Loci::multiMap ;
using Loci::const_multiMap ;
using Loci::inverseMap ;

using Loci::store ;
using Loci::const_store ;

using Loci::storeVec ;
using Loci::const_storeVec ;
using Loci::storeMat ;
using Loci::const_storeMat ;
using Loci::pivot_type ;

using Loci::param ;
using Loci::const_param ;

using Loci::constraint ;

using Loci::path ;

using Loci::unit_rule ;
using Loci::apply_rule ;
using Loci::pointwise_rule ;
using Loci::register_rule ;


using Loci::sequence ;
using Loci::do_loop ;
using Loci::Entity ;

using Loci::global_rule_list ;
using Loci::rule_db ;

using Loci::fact_db ;
using Loci::executeP ;
using Loci::create_execution_schedule ;

using Loci::set_fpe_abort ;

using Loci::options_list ;
#endif





