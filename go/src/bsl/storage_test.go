package bsl

import (
	"testing"
)

//test of storage creation
func TestStorageCreation(t *testing.T) {
	st := NewStorage(2, 5, 7)
	st.Set(1, 4, 6, complex(1, 2))
	
	if (st.At(1,4,6) != complex(1,2)) {
		t.Error("set/get storage error!")
	}
}
