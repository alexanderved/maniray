#ifndef _MR_MISC_H
#define _MR_MISC_H

#define MR_FAILURE 0
#define MR_SUCCESS 1

#define MR_UNUSED(x) ((void)(x))

#define MR_INVALID_INDEX -1

#define MR_DEFINE_CALLBACK(name, ret_type, ...) \
    typedef ret_type (*name ## _fn)(__VA_ARGS__, void *userdata); \
    typedef struct name ## _cb { \
        name ## _fn fn; \
        void *userdata; \
    } name ## _cb; \
    static inline name ## _cb name ## _cb_create(name ## _fn fn, void *userdata) { \
        return (name ## _cb) { fn, userdata }; \
    } \
    static inline name ## _cb name ## _cb_null() { \
        return (name ## _cb) { NULL, NULL }; \
    }

#endif // _MR_MISC_H