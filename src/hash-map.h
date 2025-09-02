#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define INITIAL_CAPACITY 256
#define LOAD_FACTOR_THRESHOLD 0.7
#define UB 31 //upperbound
#define VECTOR_COMPARE_FULL_MATCH 0x0101010101010101ULL

static const uint64_t base_to_bits[128] = { ['A'] = 0, ['C'] = 1, ['T'] = 2, ['G'] = 3 };

// -------------------- Custom Instructions --------------------

static inline uint64_t accelerator_reduce_sum(uint64_t a, uint64_t b) {
    uint64_t result;
    asm volatile (
        "ld x10, %[a0]\n"
        "ld x11, %[a1]\n"
        "accelerator_reduce_sum %[s], x10, x11\n"
        : [s]"=r"(result)
        : [a0]"m"(a), [a1]"m"(b)
        : "x10", "x11"
    );
    return result;
}

// same -> result = 0
static inline uint64_t vector_compare(uint64_t a, uint64_t b) {
    uint64_t result;
    asm volatile (
        "ld x10, %[a0]\n\t"
        "ld x11, %[a1]\n\t"
        "vector_compare %[s], x10, x11\n\t"
        : [s]"=r"(result)
        : [a0]"m"(a), [a1]"m"(b)
        : "x10", "x11"
    );
    return result;
}

// -------------------- Structs --------------------

typedef struct {
    int *data;
    int size;
    int capacity;
} IntArray;

typedef struct Entry {
    char *key;
    IntArray value;
    struct Entry *next;
} Entry;

typedef struct {
    Entry **buckets;
    int size;       // number of entries
    int capacity;   // number of buckets
} HashMap;

//mastermap************************************
typedef struct MasterEntry {
    uint64_t kmer_code;
    IntArray value;
    struct MasterEntry *next;
} MasterEntry;

typedef struct {
    MasterEntry **buckets;
    int size;       // number of entries
    int capacity;   // number of buckets
} SubMap;
//mastermap************************************

//treemap************************************
typedef struct TreeNode { //node
    IntArray value; //leaf node
    struct TreeMap *next; //middle node
} TreeNode;

typedef struct TreeEntry { //entry
    uint64_t key; //key
    TreeNode *child; //node
    struct TreeEntry *next; //collision
} TreeEntry;

typedef struct TreeMap {
    TreeEntry **buckets; //entry
    int size;
    int capacity;
} TreeMap;
//treemap-************************************

// -------------------- IntArray Functions --------------------

void init_int_array(IntArray *arr) {
    arr->size = 0;
    arr->capacity = 4;
    arr->data = malloc(arr->capacity * sizeof(int));
}

void add_to_int_array(IntArray *arr, int value) {
    if (arr->size >= arr->capacity) {
        arr->capacity += 4;
        arr->data = realloc(arr->data, arr->capacity * sizeof(int));
    }
    arr->data[arr->size++] = value;
}

void free_int_array(IntArray *arr) {
    if (arr->data != NULL) {
        free(arr->data);
        arr->data = NULL;
    }
    arr->size = 0;
    arr->capacity = 0;
}

// -------------------- Hash Functions --------------------

unsigned int hash_original(const char *str, int capacity) {
    unsigned long hash = 5381;
    int c;
    while ((c = *str++))
        hash = ((hash << 5) + hash) + c;
    return hash % capacity;
}

uint64_t splitmix64(uint64_t x) {
    // x = accelerator_reduce_sum(x, 0x9e3779b97f4a7c15);
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return x;
}

unsigned int hash_uint64(uint64_t key, int capacity) {
    return splitmix64(key) & (capacity - 1); // capacity must be power of 2
}

uint64_t kmer_encode_init(const char *kmer, int k) {
    uint64_t result = 0;
    for (int i = 0; i < k; i++) {
        result = (result << 2) | base_to_bits[kmer[i]];
    }
    return result;
}

uint64_t kmer_encode_roll(uint64_t prev, char new_char, int k) {
    uint64_t mask = ((uint64_t)1 << (2 * k)) - 1;
    prev = ((prev << 2) | base_to_bits[new_char]) & mask;
    return prev;
}

void kmer_encode_segs_init(const char *kmer, int k, uint64_t *codes) {
    int seg_count = (k + UB - 1) / UB;
    int i = 0;
    for (i = 0; i < (seg_count-1); i++) {
        codes[i] = kmer_encode_init(kmer + (i*UB), UB);
    }
    codes[seg_count-1] = kmer_encode_init(kmer + (i*UB), k-(i*UB));
}

void kmer_encode_segs_roll(uint64_t *codes, const char *kmer, int k, int seg_count) {
    int i = 0;
    for(i = 0; i < (seg_count-1); i++) {
        codes[i] = kmer_encode_roll(codes[i], kmer[(i+1)*UB-1], UB);
    }
    codes[i] = kmer_encode_roll(codes[i], kmer[k-1], k-(i*UB));
}

unsigned int hash(const char *str, int capacity) {
    int len = strlen(str)/32;
    uint64_t h1 = 0x0000000000000000;
    for (int i = 0; i < len; i++) {
    	uint64_t code = kmer_encode_init(str, 32);
	    h1 = h1 ^ splitmix64(code);
	    str += 32;
    }
    uint64_t code = kmer_encode_init(str, strlen(str));
    h1 = h1 ^ splitmix64(code);
    return h1 & (capacity-1);
}

void report_hashmap_stats(SubMap *map) {
    int used_buckets = 0;
    int max_chain = 0;
    long total_collisions = 0;

    for (int i = 0; i < map->capacity; i++) {
        MasterEntry *entry = map->buckets[i];
        if (entry) used_buckets++;

        int chain_length = 0;
        while (entry) {
            chain_length++;
            entry = entry->next;
        }
        if (chain_length > 1) {
            total_collisions += (chain_length - 1);
        }
        if (chain_length > max_chain) {
            max_chain = chain_length;
        }
    }

    printf("HashMap Capacity: %d\n", map->capacity);
    printf("Used Buckets: %d\n", used_buckets);
    printf("Load Factor: %.2f\n", (float)map->size / map->capacity);
    printf("Max Chain Length: %d\n", max_chain);
    printf("Total Collisions: %ld\n", total_collisions);
}

// -------------------- HashMap Functions --------------------

HashMap* create_hashmap() {
    HashMap *map = malloc(sizeof(HashMap));
    map->capacity = INITIAL_CAPACITY;
    map->size = 0;
    map->buckets = calloc(map->capacity, sizeof(Entry*));
    return map;
}

void free_entry_list(Entry *entry) {
    while (entry) {
        Entry *next = entry->next;
        free(entry->key);
        free(entry->value.data);
        free(entry);
        entry = next;
    }
}

void free_hashmap(HashMap *map) {
    for (int i = 0; i < map->capacity; i++) {
        free_entry_list(map->buckets[i]);
    }
    free(map->buckets);
    free(map);
}

// YU modified
void clean_hashmap(HashMap *map) {
    for (int i = 0; i < map->capacity; i++) {
        free_entry_list(map->buckets[i]);
    }
    free(map->buckets);
    map->buckets = calloc(map->capacity, sizeof(Entry*));
    map->size = 0;
}

// Rehashes all entries into a new bucket array
void resize_hashmap(HashMap *map) {
    int oldCapacity = map->capacity;
    Entry **oldBuckets = map->buckets;

    // map->capacity += INITIAL_CAPACITY;
    map->capacity *= 2;

    map->buckets = calloc(map->capacity, sizeof(Entry*));
    map->size = 0;

    for (int i = 0; i < oldCapacity; i++) {
        Entry *entry = oldBuckets[i];
        while (entry) {
            Entry *next = entry->next;

            // Reinsert key-value into new table
            unsigned int index = hash(entry->key, map->capacity);
            entry->next = map->buckets[index];
            map->buckets[index] = entry;

            map->size++;
	        // map->size = accelerator_reduce_sum(map->size, 1);

            entry = next;
        }
    }
    free(oldBuckets);
}

// Insert or update a key with a new value
void insert(HashMap *map, const char *key, int value) {
    double loadFactor = (double)(map->size + 1) / map->capacity;
    if (loadFactor > LOAD_FACTOR_THRESHOLD) {
        resize_hashmap(map);
    }

    unsigned int index = hash(key, map->capacity);
    Entry *entry = map->buckets[index];

    // while (entry) {
    //     if (strcmp(entry->key, key) == 0) {
    //         add_to_int_array(&entry->value, value);
    //         return;
    //     }
    //     entry = entry->next;
    // }
    while (entry) {
        uint64_t vector_comparison_result = 0;
        for (int t = 0; t < strlen(key); t+=8){
            uint64_t int1 = 0;
            uint64_t int2 = 0;
            for (int tt = t; tt < strlen(key) && tt < t+8; tt++) {
                int1 = (int1 << 8) | (entry->key)[tt];
                int2 = (int2 << 8) | (key)[tt];
            }
            vector_comparison_result = vector_compare(int1, int2);
            if (vector_comparison_result) break;
        }
        if (!vector_comparison_result) {
            add_to_int_array(&entry->value, value);
            return;
	    }
        entry = entry->next;
    }

    // Key not found, create new entry
    Entry *newEntry = malloc(sizeof(Entry));
    newEntry->key = strdup(key);
    init_int_array(&newEntry->value);
    add_to_int_array(&newEntry->value, value);
    newEntry->next = map->buckets[index];
    map->buckets[index] = newEntry;
    map->size++;
}

void print_values(HashMap *map, const char *key) {
    unsigned int index = hash(key, map->capacity);
    Entry *entry = map->buckets[index];

    while (entry) {
        if (strcmp(entry->key, key) == 0) {
            printf("%s : ", key);
            for (int i = 0; i < entry->value.size; i++) {
                printf("%d ", entry->value.data[i]);
            }
            printf("\n");
            return;
        }
        entry = entry->next;
    }

    // printf("Key \"%s\" not found.\n", key);
}

void get_values(HashMap *map, const char *key, IntArray *values) {
    unsigned int index = hash(key, map->capacity);
    Entry *entry = map->buckets[index];

    while (entry) {
        if (strcmp(entry->key, key) == 0) {
            for (int i = 0; i < entry->value.size; i++) {
                add_to_int_array(values, entry->value.data[i]);
            }
        }
        entry = entry->next;
    }
}

void write_hashmap_to_file(HashMap *map, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Could not open file for writing");
        return;
    }

    for (int i = 0; i < map->capacity; i++) {
        Entry *entry = map->buckets[i];
        while (entry) {
            fprintf(file, "%s :", entry->key);
            for (int j = 0; j < entry->value.size; j++) {
                fprintf(file, " %d", entry->value.data[j]);
            }
            fprintf(file, "\n");
            entry = entry->next;
        }
    }
    fclose(file);
}

void read_hashmap_from_file(HashMap *map, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file for reading");
        return;
    }

    char line[10240];

    int i;
    while (fgets(line, sizeof(line), file)) {
        char *key = strtok(line, " :\t\n");
        if (!key) continue;

        char *token = strtok(NULL, " \t\n");
	    i = 0;
        while (token) {
            int value = atoi(token);
	        if (i > 0)
                insert(map, key, value);  // Reuse existing insert logic
            token = strtok(NULL, " \t\n");
	        ++i;
        }
    }

    fclose(file);
}

// -------------------- MasterMap Functions --------------------

SubMap* create_submap() {
    SubMap *map = malloc(sizeof(SubMap));
    map->capacity = INITIAL_CAPACITY;
    map->size = 0;
    map->buckets = calloc(map->capacity, sizeof(MasterEntry*));
    return map;
}

void free_Masterentry_list(MasterEntry *entry) {
    while (entry) {
        MasterEntry *next = entry->next;
        free(entry->value.data);
        free(entry);
        entry = next;
    }
}

void free_Subhashmap(SubMap *map) {
    for (int i = 0; i < map->capacity; i++) {
        free_Masterentry_list(map->buckets[i]);
    }
    free(map->buckets);
    free(map);
}

// Rehashes all entries into a new bucket array
void resize_Submap(SubMap *map) {
    int oldCapacity = map->capacity;
    MasterEntry **oldBuckets = map->buckets;

    // map->capacity += INITIAL_CAPACITY;
    map->capacity *= 2;
    map->buckets = calloc(map->capacity, sizeof(MasterEntry*));
    map->size = 0;

    for (int i = 0; i < oldCapacity; i++) {
        MasterEntry *entry = oldBuckets[i];
        while (entry) {
            MasterEntry *next = entry->next;

            // Reinsert key-value into new table
            unsigned int index = hash_uint64(entry->kmer_code, map->capacity);
            entry->next = map->buckets[index];
            map->buckets[index] = entry;

            map->size++;
            entry = next;
        }
    }
    free(oldBuckets);
}

// Insert or update a key with a new value
void insert_Submap(SubMap *map, const uint64_t kmer_code, int value) {
    double loadFactor = (double)(map->size + 1) / map->capacity;
    if (loadFactor > LOAD_FACTOR_THRESHOLD) {
        resize_Submap(map);
    }
    unsigned int index = hash_uint64(kmer_code, map->capacity);
    MasterEntry *entry = map->buckets[index];

    while (entry) {
        if (entry->kmer_code == kmer_code) {
            add_to_int_array(&entry->value, value);
            return;
        }
        entry = entry->next;
    }

    // Key not found, create new entry
    MasterEntry *newEntry = malloc(sizeof(MasterEntry));
    newEntry->kmer_code = kmer_code;
    init_int_array(&newEntry->value);
    add_to_int_array(&newEntry->value, value);
    newEntry->next = map->buckets[index];
    map->buckets[index] = newEntry;
    map->size++;
}

void get_values_Submap(SubMap *map, const uint64_t kmer_code, int *matched_transcripts_counts, int *max_val, int *sec_val) {
    unsigned int index = hash_uint64(kmer_code, map->capacity);
    MasterEntry *entry = map->buckets[index];

    while (entry) {
        if (!vector_compare(entry->kmer_code, kmer_code)) {
            for (int i = 0; i < entry->value.size; i++) {
                int idx = entry->value.data[i];
                if(matched_transcripts_counts[idx] < (*max_val)) {
                    if(matched_transcripts_counts[idx] == (*sec_val)) {
                        (*sec_val)++;
                    }
                }
                else if(matched_transcripts_counts[idx] == (*max_val)) {
                    (*max_val)++;
                }
                matched_transcripts_counts[idx]++;
            }
            break;
        }
        entry = entry->next;
    }
}

void read_Submap_from_file(SubMap *map, const char *filename, int k) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file for reading");
        return;
    }

    char line[10240];
    int i;
    while (fgets(line, sizeof(line), file)) {

        char *key = strtok(line, " :\t\n");
        if (!key) continue;
        uint64_t encode_key = kmer_encode_init(key, k);

        char *token = strtok(NULL, " \t\n");
        i = 0;
        while (token) {
            int value = atoi(token);
            if (i > 0)
                insert_Submap(map, encode_key, value);
            token = strtok(NULL, " \t\n");
            ++i;
        }
    }
    fclose(file);
}

// -------------------- TreeMap Functions --------------------

TreeMap* create_treemap() {
    TreeMap *map = malloc(sizeof(TreeMap));
    map->capacity = INITIAL_CAPACITY;
    map->size = 0;
    map->buckets = calloc(map->capacity, sizeof(TreeEntry*));
    return map;
}

void free_treemap(TreeMap *map);

void free_Treeentry_list(TreeEntry *entry) {
    while (entry) {
        TreeEntry *next = entry->next;
        if(entry->child) {
            free(entry->child->value.data);
            if(entry->child->next) {
                free_treemap(entry->child->next);
            }
            free(entry->child);
        }
        free(entry);
        entry = next;
    }
}

void free_treemap(TreeMap *map) {
    for (int i = 0; i < map->capacity; i++) {
        free_Treeentry_list(map->buckets[i]);
    }
    free(map->buckets);
    free(map);
}

void resize_treemap(TreeMap *map) {
    int old_capacity = map->capacity;
    TreeEntry **old_buckets = map->buckets;

    map->capacity *= 2;
    map->buckets = calloc(map->capacity, sizeof(TreeEntry*));
    map->size = 0;

    for (int i = 0; i < old_capacity; i++) {
        TreeEntry *entry = old_buckets[i];
        while (entry) {
            TreeEntry *next = entry->next;

            unsigned int index = hash_uint64(entry->key, map->capacity);

            entry->next = map->buckets[index];
            map->buckets[index] = entry;

            map->size++;
            entry = next;
        }
    }
    free(old_buckets);
}

void insert_Treemap(TreeMap *map, const uint64_t *codes, int seg_count, int value) {
    double load_factor = (double)(map->size + 1) / map->capacity;
    if (load_factor > LOAD_FACTOR_THRESHOLD) {
        resize_treemap(map);
    }

    TreeMap *current_map = map;

    for (int i = 0; i < seg_count; ++i) {
        uint64_t code = codes[i];
        unsigned int index = hash_uint64(code, current_map->capacity);
        TreeEntry *entry = current_map->buckets[index];

        while (entry) {
            if (entry->key == code)
                break;
            entry = entry->next;
        }

        if (!entry) {
            entry = malloc(sizeof(TreeEntry));
            entry->key = code;
            entry->next = current_map->buckets[index];
            current_map->buckets[index] = entry;
            current_map->size++;

            entry->child = malloc(sizeof(TreeNode));
            entry->child->next = NULL;
            init_int_array(&entry->child->value);
        }

        if (i == seg_count - 1) { // leaf node
            add_to_int_array(&entry->child->value, value);
        }
        else {
            if (!entry->child->next) {
                entry->child->next = create_treemap();
            }
            current_map = entry->child->next;
        }
    }
}

void get_values_Treemap(TreeMap *map, const uint64_t *codes, int seg_count, int *matched_transcripts_counts, int *max_val, int *sec_val) {
    TreeMap *current_map = map;
    TreeEntry *entry = NULL;

    for (int i = 0; i < seg_count; ++i) {
        uint64_t code = codes[i];
        unsigned int index = hash_uint64(code, current_map->capacity);
        entry = current_map->buckets[index];

        while (entry) {
            if (!vector_compare(entry->key, code))
                break;
            entry = entry->next;
        }

        if (!entry) {
            return;  // key not found, nothing to update
        }

        if (i < seg_count - 1) {
            if (!entry->child->next) return; // middle node but no next
            current_map = entry->child->next;
        }
    }

    // in the last layer, find leaf node
    if (entry) {
        IntArray *arr = &entry->child->value;
        for (int i = 0; i < arr->size; i++) {
            int idx = arr->data[i];
            if (matched_transcripts_counts[idx] < (*max_val)) {
                if (matched_transcripts_counts[idx] == (*sec_val)) {
                    (*sec_val)++;
                }
            } else if (matched_transcripts_counts[idx] == (*max_val)) {
                (*max_val)++;
            }
            
            matched_transcripts_counts[idx]++;
        }
    }
}

void read_Treemap_from_file(TreeMap *map, const char *filename, int k) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file for reading");
        return;
    }

    int seg_count = (k + UB - 1) / UB; // ceiling
    uint64_t *codes = malloc(sizeof(uint64_t) * seg_count);

    char line[10240];
    int i;

    while (fgets(line, sizeof(line), file)) {
        char *key = strtok(line, " :\t\n");
        if (!key) continue;

        kmer_encode_segs_init(key, k, codes);

        char *token = strtok(NULL, " \t\n");
        i = 0;
        while (token) {
            int value = atoi(token);
            if (i > 0)
                insert_Treemap(map, codes, seg_count, value);
            token = strtok(NULL, " \t\n");
            ++i;
        }
    }
    free(codes);
    fclose(file);
}

void print_values_Treemap(TreeMap *map, const char *kmer, int k) {
    int seg_count = (k + UB - 1) / UB;
    uint64_t *codes = malloc(sizeof(uint64_t) * seg_count);
    kmer_encode_segs_init(kmer, k, codes);

    TreeMap *current_map = map;
    TreeEntry *entry = NULL;

    for (int i = 0; i < seg_count; ++i) {
        uint64_t code = codes[i];
        unsigned int index = hash_uint64(code, current_map->capacity);
        entry = current_map->buckets[index];

        while (entry) {
            if (entry->key == code)
                break;
            entry = entry->next;
        }

        if (!entry) {
            // printf("Key \"%s\" not found.\n", kmer);
            free(codes);
            return;
        }

        if (i < seg_count - 1) {
            if (!entry->child->next) {
                // printf("Key \"%s\" not found (incomplete tree).\n", kmer);
                free(codes);
                return;
            }
            current_map = entry->child->next;
        }
    }

    if (entry) {
        IntArray *arr = &entry->child->value;
        for(int i=0; i<k; i++) {
            printf("%c", kmer[i]);
        }
        printf(" : ");
        for (int i = 0; i < arr->size; i++) {
            printf("%d ", arr->data[i]);
        }
        printf("\n");
    }

    free(codes);
}
