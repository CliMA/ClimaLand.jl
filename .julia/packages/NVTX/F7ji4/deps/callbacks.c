#include <stdint.h>
#include <nvtx3/nvToolsExt.h>
#include "callbacks.h"

nvtxDomainHandle_t julia_domain = 0;
nvtxStringHandle_t gc_message = 0;
nvtxStringHandle_t gc_alloc_message = 0;
nvtxStringHandle_t gc_free_message = 0;
uint32_t gc_color = 0;
uint32_t gc_alloc_color = 0;
uint32_t gc_free_color = 0;

// https://github.com/JuliaLang/julia/blob/v1.8.3/src/julia_gcext.h#L20-L22
extern void nvtx_julia_gc_cb_pre(int full) {
  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = gc_color;
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_REGISTERED;
  eventAttrib.message.registered = gc_message;
  eventAttrib.category = (uint32_t) (1 + full);
  nvtxDomainRangePushEx(julia_domain, &eventAttrib);
}

extern void nvtx_julia_gc_cb_post(int full) {
  nvtxDomainRangePop(julia_domain);
}

// https://github.com/JuliaLang/julia/blob/v1.8.3/src/julia_gcext.h#L24-L26
extern void nvtx_julia_gc_cb_alloc(void *ptr, size_t size) {
  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = gc_alloc_color;
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_REGISTERED;
  eventAttrib.message.registered = gc_alloc_message;
  eventAttrib.payloadType = NVTX_PAYLOAD_TYPE_UNSIGNED_INT64;
  eventAttrib.payload.ullValue = (uint64_t) size;
  nvtxDomainMarkEx(julia_domain, &eventAttrib);
}
extern void nvtx_julia_gc_cb_free(void *ptr) {
  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = gc_free_color;
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_REGISTERED;
  eventAttrib.message.registered = gc_free_message;
  nvtxDomainMarkEx(julia_domain, &eventAttrib);
}
