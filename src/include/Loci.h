#ifndef LOCI_H
#define LOCI_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Loci_version.h>
#include <Loci_types.h>
#include <data_traits.h>
#include <Map.h>
#include <MapVec.h>
#include <multiMap.h>
#include <store.h>
#include <storeVec.h>
#include <parameter.h>
#include <constraint.h>
#include <rule.h>
#include <Tools/options_list.h>
#include <Tools/fpe.h>
#include <scheduler.h>
#include <distribute.h>
#include <multiStore.h>
#include <storeMat.h>
#include <mod_db.h>

using Loci::entitySet ;
using Loci::create_entitySet ;
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
using Loci::multiStore ;
using Loci::const_multiStore ;
using Loci::storeVec ;
using Loci::const_storeVec ;
using Loci::Vect ;
using Loci::Mat ;
using Loci::const_Vect ;
using Loci::const_Mat ;
using Loci::Scalar ;
using Loci::mk_Scalar ;
using Loci::pivot_type ;

using Loci::param ;
using Loci::const_param ;

using Loci::constraint ;

using Loci::unit_rule ;
using Loci::apply_rule ;
using Loci::pointwise_rule ;
using Loci::singleton_rule ;
using Loci::register_rule ;


using Loci::storeMat ;
using Loci::const_storeMat ;

#include <DStore.h>
using Loci::dstore ;
using Loci::const_dstore ;

#include <DStoreVec.h>
using Loci::dstoreVec;
using Loci::const_dstoreVec;

//#include <DStoreMat.h>
//using Loci::dstoreMat;
//using Loci::const_dstoreMat;


#include <DMultiStore.h>
using Loci::dmultiStore;
using Loci::const_dmultiStore;

#include <DMultiMap.h>
using Loci:: dmultiMap;
using Loci:: const_dmultiMap;

#include <DMap.h>
using Loci:: dMap;
using Loci:: const_dMap;

#include <DMapVec.h>
using Loci:: dMapVec;
using Loci:: const_dMapVec;

using Loci::sequence ;
using Loci::create_sequence ;
using Loci::do_loop ;
using Loci::Entity ;

using Loci::global_rule_list ;
using Loci::rule_db ;

using Loci::fact_db ;
using Loci::executeP ;
using Loci::create_execution_schedule ;

using Loci::set_fpe_abort ;

using Loci::options_list ;

using Loci::data_schema_traits ;
using Loci::data_schema_converter_traits ;

using Loci::vector3d ;
using Loci::vector2d ;
using Loci::norm ;
using Loci::dot ;
using Loci::cross ;

using Loci::Array ;
using Loci::mod ;
using Loci::load_module ;
using Loci::rule_impl_list ;
#endif





