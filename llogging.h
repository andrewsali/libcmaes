
#ifndef LLOGGING_H
#define LLOGGING_H

#include "libcmaes_config.h"

#ifdef HAVE_GLOG // HAVE_LIB_GLOG
#include <glog/logging.h>
#else
#include <iostream>

namespace libcmaes
{
#define INFO "INFO"
#define WARNING "WARNING"
#define ERROR "ERROR"
#define FATAL "FATAL"

static std::ostream nullstream(0);

inline std::ostream& LOG(const std::string &severity,std::ostream &out=std::cout)
{
  out << severity << " - ";
  return out;
}

inline std::ostream& LOG_IF(const std::string &severity,const bool &condition,std::ostream &out=std::cout)
{
  if (condition)
    return LOG(severity,out);
  else return nullstream;
}
}
#endif
#endif
