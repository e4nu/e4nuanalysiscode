// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME _ROOT_DICT_E4NuPhysics

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
#include "EventHolderI.h"

// Header files passed via #pragma extra_include

namespace e4nu {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *e4nu_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("e4nu", 0 /*version*/, "EventHolderI.h", 16,
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

namespace ROOT {
   static TClass *e4nucLcLEventHolderI_Dictionary();
   static void e4nucLcLEventHolderI_TClassManip(TClass*);
   static void *new_e4nucLcLEventHolderI(void *p = 0);
   static void *newArray_e4nucLcLEventHolderI(Long_t size, void *p);
   static void delete_e4nucLcLEventHolderI(void *p);
   static void deleteArray_e4nucLcLEventHolderI(void *p);
   static void destruct_e4nucLcLEventHolderI(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::e4nu::EventHolderI*)
   {
      ::e4nu::EventHolderI *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::e4nu::EventHolderI));
      static ::ROOT::TGenericClassInfo 
         instance("e4nu::EventHolderI", "EventHolderI.h", 17,
                  typeid(::e4nu::EventHolderI), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &e4nucLcLEventHolderI_Dictionary, isa_proxy, 0,
                  sizeof(::e4nu::EventHolderI) );
      instance.SetNew(&new_e4nucLcLEventHolderI);
      instance.SetNewArray(&newArray_e4nucLcLEventHolderI);
      instance.SetDelete(&delete_e4nucLcLEventHolderI);
      instance.SetDeleteArray(&deleteArray_e4nucLcLEventHolderI);
      instance.SetDestructor(&destruct_e4nucLcLEventHolderI);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::e4nu::EventHolderI*)
   {
      return GenerateInitInstanceLocal((::e4nu::EventHolderI*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::e4nu::EventHolderI*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *e4nucLcLEventHolderI_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::e4nu::EventHolderI*)0x0)->GetClass();
      e4nucLcLEventHolderI_TClassManip(theClass);
   return theClass;
   }

   static void e4nucLcLEventHolderI_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_e4nucLcLEventHolderI(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::e4nu::EventHolderI : new ::e4nu::EventHolderI;
   }
   static void *newArray_e4nucLcLEventHolderI(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::e4nu::EventHolderI[nElements] : new ::e4nu::EventHolderI[nElements];
   }
   // Wrapper around operator delete
   static void delete_e4nucLcLEventHolderI(void *p) {
      delete ((::e4nu::EventHolderI*)p);
   }
   static void deleteArray_e4nucLcLEventHolderI(void *p) {
      delete [] ((::e4nu::EventHolderI*)p);
   }
   static void destruct_e4nucLcLEventHolderI(void *p) {
      typedef ::e4nu::EventHolderI current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::e4nu::EventHolderI

namespace {
  void TriggerDictionaryInitialization_libE4NuPhysics_Impl() {
    static const char* headers[] = {
"EventHolderI.h",
0
    };
    static const char* includePaths[] = {
"/grid/fermiapp/products/larsoft/root/v6_12_06a/Linux64bit+3.10-2.17-e17-debug/include",
"/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/src/",
"../include",
"/grid/fermiapp/products/larsoft/root/v6_12_06a/Linux64bit+3.10-2.17-e17-debug/include",
"/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/src/",
"/grid/fermiapp/products/larsoft/root/v6_12_06a/Linux64bit+3.10-2.17-e17-debug/include",
"/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/src/physics/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libE4NuPhysics dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace e4nu{class __attribute__((annotate("$clingAutoload$EventHolderI.h")))  EventHolderI;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libE4NuPhysics dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "EventHolderI.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"e4nu::EventHolderI", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libE4NuPhysics",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libE4NuPhysics_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libE4NuPhysics_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libE4NuPhysics() {
  TriggerDictionaryInitialization_libE4NuPhysics_Impl();
}
