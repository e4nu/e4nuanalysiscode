// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME _ROOT_DICT_E4NuConf

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "FiducialCutI.h"
#include "AccpetanceMapsI.h"
#include "AnalysisCutsI.h"
#include "ConfigurablesI.h"

// Header files passed via #pragma extra_include

namespace e4nu {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *e4nu_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("e4nu", 0 /*version*/, "FiducialCutI.h", 8,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &e4nu_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *e4nu_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace e4nu {
   namespace conf {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *e4nucLcLconf_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("e4nu::conf", 0 /*version*/, "FiducialCutI.h", 9,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &e4nucLcLconf_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *e4nucLcLconf_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}
}

namespace {
  void TriggerDictionaryInitialization_libE4NuConf_Impl() {
    static const char* headers[] = {
"FiducialCutI.h",
"AccpetanceMapsI.h",
"AnalysisCutsI.h",
"ConfigurablesI.h",
0
    };
    static const char* includePaths[] = {
"/grid/fermiapp/products/larsoft/root/v6_12_06a/Linux64bit+3.10-2.17-e17-debug/include",
"/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/src/",
"../include",
"/grid/fermiapp/products/larsoft/root/v6_12_06a/Linux64bit+3.10-2.17-e17-debug/include",
"/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/src/",
"/grid/fermiapp/products/larsoft/root/v6_12_06a/Linux64bit+3.10-2.17-e17-debug/include",
"/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/src/conf/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libE4NuConf dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libE4NuConf dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "FiducialCutI.h"
#include "AccpetanceMapsI.h"
#include "AnalysisCutsI.h"
#include "ConfigurablesI.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libE4NuConf",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libE4NuConf_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libE4NuConf_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libE4NuConf() {
  TriggerDictionaryInitialization_libE4NuConf_Impl();
}
