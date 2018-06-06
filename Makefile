CC = gcc
CFLAGS = -O3 -Wno-write-strings -fPIC $(shell pkg-config --cflags glib-2.0)
LDFLAGS = -shared $(shell pkg-config --libs glib-2.0)

TARGET  = src/PPCAS.so
SOURCES = $(shell echo src/*.c)
OBJECTS = $(SOURCES:.c=.o)

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)
	
$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)