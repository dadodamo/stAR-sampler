// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: paramdata.proto

#include "paramdata.pb.h"

#include <algorithm>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/wire_format_lite.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
extern PROTOBUF_INTERNAL_EXPORT_paramdata_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_matrix_paramdata_2eproto;
extern PROTOBUF_INTERNAL_EXPORT_paramdata_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_vector_paramdata_2eproto;
namespace sampler_data {
class vectorDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<vector> _instance;
} _vector_default_instance_;
class matrixDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<matrix> _instance;
} _matrix_default_instance_;
class samplesDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<samples> _instance;
} _samples_default_instance_;
}  // namespace sampler_data
static void InitDefaultsscc_info_matrix_paramdata_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::sampler_data::_matrix_default_instance_;
    new (ptr) ::sampler_data::matrix();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_matrix_paramdata_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_matrix_paramdata_2eproto}, {
      &scc_info_vector_paramdata_2eproto.base,}};

static void InitDefaultsscc_info_samples_paramdata_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::sampler_data::_samples_default_instance_;
    new (ptr) ::sampler_data::samples();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<2> scc_info_samples_paramdata_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 2, 0, InitDefaultsscc_info_samples_paramdata_2eproto}, {
      &scc_info_matrix_paramdata_2eproto.base,
      &scc_info_vector_paramdata_2eproto.base,}};

static void InitDefaultsscc_info_vector_paramdata_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::sampler_data::_vector_default_instance_;
    new (ptr) ::sampler_data::vector();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_vector_paramdata_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_vector_paramdata_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_paramdata_2eproto[3];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_paramdata_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_paramdata_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_paramdata_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::sampler_data::vector, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::sampler_data::vector, vec_value_),
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::sampler_data::matrix, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::sampler_data::matrix, vec_),
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, o_),
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, beta_),
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, mu_0_),
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, rho_),
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, sigma_w_),
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, sigma_0_),
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, sigma_eps_),
  PROTOBUF_FIELD_OFFSET(::sampler_data::samples, phi_),
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, -1, sizeof(::sampler_data::vector)},
  { 6, -1, sizeof(::sampler_data::matrix)},
  { 12, -1, sizeof(::sampler_data::samples)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::sampler_data::_vector_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::sampler_data::_matrix_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::sampler_data::_samples_default_instance_),
};

const char descriptor_table_protodef_paramdata_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\017paramdata.proto\022\014sampler_data\"\033\n\006vecto"
  "r\022\021\n\tvec_value\030\001 \003(\001\"+\n\006matrix\022!\n\003vec\030\001 "
  "\003(\0132\024.sampler_data.vector\"\301\001\n\007samples\022\037\n"
  "\001o\030\001 \003(\0132\024.sampler_data.matrix\022\"\n\004beta\030\002"
  " \003(\0132\024.sampler_data.vector\022\"\n\004mu_0\030\003 \003(\013"
  "2\024.sampler_data.vector\022\013\n\003rho\030\004 \003(\001\022\017\n\007s"
  "igma_w\030\005 \003(\001\022\017\n\007sigma_0\030\006 \003(\001\022\021\n\tsigma_e"
  "ps\030\007 \003(\001\022\013\n\003phi\030\010 \003(\001"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_paramdata_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_paramdata_2eproto_sccs[3] = {
  &scc_info_matrix_paramdata_2eproto.base,
  &scc_info_samples_paramdata_2eproto.base,
  &scc_info_vector_paramdata_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_paramdata_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_paramdata_2eproto = {
  false, false, descriptor_table_protodef_paramdata_2eproto, "paramdata.proto", 301,
  &descriptor_table_paramdata_2eproto_once, descriptor_table_paramdata_2eproto_sccs, descriptor_table_paramdata_2eproto_deps, 3, 0,
  schemas, file_default_instances, TableStruct_paramdata_2eproto::offsets,
  file_level_metadata_paramdata_2eproto, 3, file_level_enum_descriptors_paramdata_2eproto, file_level_service_descriptors_paramdata_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_paramdata_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_paramdata_2eproto)), true);
namespace sampler_data {

// ===================================================================

class vector::_Internal {
 public:
};

vector::vector(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  vec_value_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:sampler_data.vector)
}
vector::vector(const vector& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      vec_value_(from.vec_value_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:sampler_data.vector)
}

void vector::SharedCtor() {
}

vector::~vector() {
  // @@protoc_insertion_point(destructor:sampler_data.vector)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void vector::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void vector::ArenaDtor(void* object) {
  vector* _this = reinterpret_cast< vector* >(object);
  (void)_this;
}
void vector::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void vector::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const vector& vector::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_vector_paramdata_2eproto.base);
  return *internal_default_instance();
}


void vector::Clear() {
// @@protoc_insertion_point(message_clear_start:sampler_data.vector)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  vec_value_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* vector::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated double vec_value = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 9)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_vec_value(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
            ptr += sizeof(double);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<9>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_vec_value(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* vector::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:sampler_data.vector)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated double vec_value = 1;
  for (int i = 0, n = this->_internal_vec_value_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(1, this->_internal_vec_value(i), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:sampler_data.vector)
  return target;
}

size_t vector::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:sampler_data.vector)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated double vec_value = 1;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_vec_value_size());
    size_t data_size = 8UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_vec_value_size());
    total_size += data_size;
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void vector::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:sampler_data.vector)
  GOOGLE_DCHECK_NE(&from, this);
  const vector* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<vector>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:sampler_data.vector)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:sampler_data.vector)
    MergeFrom(*source);
  }
}

void vector::MergeFrom(const vector& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:sampler_data.vector)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  vec_value_.MergeFrom(from.vec_value_);
}

void vector::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:sampler_data.vector)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void vector::CopyFrom(const vector& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:sampler_data.vector)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool vector::IsInitialized() const {
  return true;
}

void vector::InternalSwap(vector* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  vec_value_.InternalSwap(&other->vec_value_);
}

::PROTOBUF_NAMESPACE_ID::Metadata vector::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

class matrix::_Internal {
 public:
};

matrix::matrix(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  vec_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:sampler_data.matrix)
}
matrix::matrix(const matrix& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      vec_(from.vec_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:sampler_data.matrix)
}

void matrix::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_matrix_paramdata_2eproto.base);
}

matrix::~matrix() {
  // @@protoc_insertion_point(destructor:sampler_data.matrix)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void matrix::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void matrix::ArenaDtor(void* object) {
  matrix* _this = reinterpret_cast< matrix* >(object);
  (void)_this;
}
void matrix::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void matrix::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const matrix& matrix::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_matrix_paramdata_2eproto.base);
  return *internal_default_instance();
}


void matrix::Clear() {
// @@protoc_insertion_point(message_clear_start:sampler_data.matrix)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  vec_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* matrix::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .sampler_data.vector vec = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_vec(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<10>(ptr));
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* matrix::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:sampler_data.matrix)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .sampler_data.vector vec = 1;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_vec_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, this->_internal_vec(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:sampler_data.matrix)
  return target;
}

size_t matrix::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:sampler_data.matrix)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .sampler_data.vector vec = 1;
  total_size += 1UL * this->_internal_vec_size();
  for (const auto& msg : this->vec_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void matrix::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:sampler_data.matrix)
  GOOGLE_DCHECK_NE(&from, this);
  const matrix* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<matrix>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:sampler_data.matrix)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:sampler_data.matrix)
    MergeFrom(*source);
  }
}

void matrix::MergeFrom(const matrix& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:sampler_data.matrix)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  vec_.MergeFrom(from.vec_);
}

void matrix::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:sampler_data.matrix)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void matrix::CopyFrom(const matrix& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:sampler_data.matrix)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool matrix::IsInitialized() const {
  return true;
}

void matrix::InternalSwap(matrix* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  vec_.InternalSwap(&other->vec_);
}

::PROTOBUF_NAMESPACE_ID::Metadata matrix::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

class samples::_Internal {
 public:
};

samples::samples(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  o_(arena),
  beta_(arena),
  mu_0_(arena),
  rho_(arena),
  sigma_w_(arena),
  sigma_0_(arena),
  sigma_eps_(arena),
  phi_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:sampler_data.samples)
}
samples::samples(const samples& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      o_(from.o_),
      beta_(from.beta_),
      mu_0_(from.mu_0_),
      rho_(from.rho_),
      sigma_w_(from.sigma_w_),
      sigma_0_(from.sigma_0_),
      sigma_eps_(from.sigma_eps_),
      phi_(from.phi_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:sampler_data.samples)
}

void samples::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_samples_paramdata_2eproto.base);
}

samples::~samples() {
  // @@protoc_insertion_point(destructor:sampler_data.samples)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void samples::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void samples::ArenaDtor(void* object) {
  samples* _this = reinterpret_cast< samples* >(object);
  (void)_this;
}
void samples::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void samples::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const samples& samples::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_samples_paramdata_2eproto.base);
  return *internal_default_instance();
}


void samples::Clear() {
// @@protoc_insertion_point(message_clear_start:sampler_data.samples)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  o_.Clear();
  beta_.Clear();
  mu_0_.Clear();
  rho_.Clear();
  sigma_w_.Clear();
  sigma_0_.Clear();
  sigma_eps_.Clear();
  phi_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* samples::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .sampler_data.matrix o = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_o(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<10>(ptr));
        } else goto handle_unusual;
        continue;
      // repeated .sampler_data.vector beta = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 18)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_beta(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<18>(ptr));
        } else goto handle_unusual;
        continue;
      // repeated .sampler_data.vector mu_0 = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 26)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_mu_0(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<26>(ptr));
        } else goto handle_unusual;
        continue;
      // repeated double rho = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 33)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_rho(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
            ptr += sizeof(double);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<33>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 34) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_rho(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated double sigma_w = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 41)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_sigma_w(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
            ptr += sizeof(double);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<41>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 42) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_sigma_w(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated double sigma_0 = 6;
      case 6:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 49)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_sigma_0(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
            ptr += sizeof(double);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<49>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 50) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_sigma_0(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated double sigma_eps = 7;
      case 7:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 57)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_sigma_eps(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
            ptr += sizeof(double);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<57>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 58) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_sigma_eps(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated double phi = 8;
      case 8:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 65)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_phi(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
            ptr += sizeof(double);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<65>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 66) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_phi(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* samples::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:sampler_data.samples)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .sampler_data.matrix o = 1;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_o_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, this->_internal_o(i), target, stream);
  }

  // repeated .sampler_data.vector beta = 2;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_beta_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, this->_internal_beta(i), target, stream);
  }

  // repeated .sampler_data.vector mu_0 = 3;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_mu_0_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(3, this->_internal_mu_0(i), target, stream);
  }

  // repeated double rho = 4;
  for (int i = 0, n = this->_internal_rho_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(4, this->_internal_rho(i), target);
  }

  // repeated double sigma_w = 5;
  for (int i = 0, n = this->_internal_sigma_w_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(5, this->_internal_sigma_w(i), target);
  }

  // repeated double sigma_0 = 6;
  for (int i = 0, n = this->_internal_sigma_0_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(6, this->_internal_sigma_0(i), target);
  }

  // repeated double sigma_eps = 7;
  for (int i = 0, n = this->_internal_sigma_eps_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(7, this->_internal_sigma_eps(i), target);
  }

  // repeated double phi = 8;
  for (int i = 0, n = this->_internal_phi_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(8, this->_internal_phi(i), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:sampler_data.samples)
  return target;
}

size_t samples::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:sampler_data.samples)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .sampler_data.matrix o = 1;
  total_size += 1UL * this->_internal_o_size();
  for (const auto& msg : this->o_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  // repeated .sampler_data.vector beta = 2;
  total_size += 1UL * this->_internal_beta_size();
  for (const auto& msg : this->beta_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  // repeated .sampler_data.vector mu_0 = 3;
  total_size += 1UL * this->_internal_mu_0_size();
  for (const auto& msg : this->mu_0_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  // repeated double rho = 4;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_rho_size());
    size_t data_size = 8UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_rho_size());
    total_size += data_size;
  }

  // repeated double sigma_w = 5;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_sigma_w_size());
    size_t data_size = 8UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_sigma_w_size());
    total_size += data_size;
  }

  // repeated double sigma_0 = 6;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_sigma_0_size());
    size_t data_size = 8UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_sigma_0_size());
    total_size += data_size;
  }

  // repeated double sigma_eps = 7;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_sigma_eps_size());
    size_t data_size = 8UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_sigma_eps_size());
    total_size += data_size;
  }

  // repeated double phi = 8;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_phi_size());
    size_t data_size = 8UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_phi_size());
    total_size += data_size;
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void samples::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:sampler_data.samples)
  GOOGLE_DCHECK_NE(&from, this);
  const samples* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<samples>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:sampler_data.samples)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:sampler_data.samples)
    MergeFrom(*source);
  }
}

void samples::MergeFrom(const samples& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:sampler_data.samples)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  o_.MergeFrom(from.o_);
  beta_.MergeFrom(from.beta_);
  mu_0_.MergeFrom(from.mu_0_);
  rho_.MergeFrom(from.rho_);
  sigma_w_.MergeFrom(from.sigma_w_);
  sigma_0_.MergeFrom(from.sigma_0_);
  sigma_eps_.MergeFrom(from.sigma_eps_);
  phi_.MergeFrom(from.phi_);
}

void samples::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:sampler_data.samples)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void samples::CopyFrom(const samples& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:sampler_data.samples)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool samples::IsInitialized() const {
  return true;
}

void samples::InternalSwap(samples* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  o_.InternalSwap(&other->o_);
  beta_.InternalSwap(&other->beta_);
  mu_0_.InternalSwap(&other->mu_0_);
  rho_.InternalSwap(&other->rho_);
  sigma_w_.InternalSwap(&other->sigma_w_);
  sigma_0_.InternalSwap(&other->sigma_0_);
  sigma_eps_.InternalSwap(&other->sigma_eps_);
  phi_.InternalSwap(&other->phi_);
}

::PROTOBUF_NAMESPACE_ID::Metadata samples::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace sampler_data
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::sampler_data::vector* Arena::CreateMaybeMessage< ::sampler_data::vector >(Arena* arena) {
  return Arena::CreateMessageInternal< ::sampler_data::vector >(arena);
}
template<> PROTOBUF_NOINLINE ::sampler_data::matrix* Arena::CreateMaybeMessage< ::sampler_data::matrix >(Arena* arena) {
  return Arena::CreateMessageInternal< ::sampler_data::matrix >(arena);
}
template<> PROTOBUF_NOINLINE ::sampler_data::samples* Arena::CreateMaybeMessage< ::sampler_data::samples >(Arena* arena) {
  return Arena::CreateMessageInternal< ::sampler_data::samples >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
