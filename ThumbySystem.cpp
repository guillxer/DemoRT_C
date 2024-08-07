#include "Arduino.h"
#include "ThumbySystem.h"
#include "pico/stdlib.h"
#include "hardware/pwm.h"
#include "SPI.h"


// Setup pins for link, buttons, and audio
void Thumby::begin(){
  // Setup link link pins
  pinMode(THUMBY_LINK_TX_PIN, OUTPUT);
  pinMode(THUMBY_LINK_RX_PIN, INPUT);
  pinMode(THUMBY_LINK_PU_PIN, OUTPUT);
  digitalWrite(THUMBY_LINK_PU_PIN, HIGH);
  
  // Setup button pins
  pinMode(THUMBY_BTN_LDPAD_PIN, INPUT_PULLUP);
  pinMode(THUMBY_BTN_RDPAD_PIN, INPUT_PULLUP);
  pinMode(THUMBY_BTN_UDPAD_PIN, INPUT_PULLUP);
  pinMode(THUMBY_BTN_DDPAD_PIN, INPUT_PULLUP);
  pinMode(THUMBY_BTN_B_PIN, INPUT_PULLUP);
  pinMode(THUMBY_BTN_A_PIN, INPUT_PULLUP);

  // Set audio pin as pwm
  gpio_set_function(THUMBY_AUDIO_PIN, GPIO_FUNC_PWM);

  // Init half-duplex UART for link cable
  Serial1.begin(115200);

  // Setup screen
  pinMode(THUMBY_CS_PIN, OUTPUT);
  pinMode(THUMBY_DC_PIN, OUTPUT);
  SPI.begin();

  // Reset the screen
  pinMode(THUMBY_RESET_PIN, OUTPUT);
  digitalWrite(THUMBY_RESET_PIN, LOW);
  delay(10);
  digitalWrite(THUMBY_RESET_PIN, HIGH);

  // Init the screen
  sendCommand(SSD1306_DISPLAYOFF);
  sendCommand(SSD1306_SETDISPLAYCLOCKDIV);
  #if GRAYSCALE_VALUES
  sendCommand(0xf0)
  #else
  sendCommand(0x80);
  #endif
  sendCommand(SSD1306_SETDISPLAYOFFSET);
  sendCommand(0x00);
  sendCommand(SSD1306_SETSTARTLINE | 0x00);
  sendCommand(SSD1306_DISPLAYALLON_RESUME);
  sendCommand(SSD1306_NORMALDISPLAY);
  sendCommand(SSD1306_CHARGEPUMP);
  sendCommand(0x14);
  sendCommand(SSD1306_MEMORYMODE);
  sendCommand(0x00);
  sendCommand(SSD1306_SEGREMAP|1);
  sendCommand(SSD1306_COMSCANDEC);
  sendCommand(SSD1306_SETCONTRAST);
  sendCommand(30);

  sendCommand(SSD1306_SETPRECHARGE);
  #if GRAYSCALE_VALUES
  sendCommand(0x11);
  #else
  sendCommand(0xF1);
  #endif

  sendCommand(SSD1306_SETVCOMDETECT);
  sendCommand(0x20);

  sendCommand(SSD1306_SETMULTIPLEX);
  sendCommand(40 - 1);

  sendCommand(SSD1306_SETCOMPINS);
  sendCommand(0x12);

  // Select internal 30uA Iref (max Iseg=240uA) during display on
  sendCommand(0xAD);
  sendCommand(0x30);

  sendCommand(SSD1306_DISPLAYON);

  // Setup screen buffer
  GraphicsBuffer::begin();
  GraphicsBuffer::setFont(thinPixel7_10ptFontInfo);
}


void Thumby::sendCommand(uint8_t command){
  digitalWrite(THUMBY_CS_PIN, 1);
  digitalWrite(THUMBY_DC_PIN, 0);
  digitalWrite(THUMBY_CS_PIN, 0);
  SPI.transfer(command);
  digitalWrite(THUMBY_CS_PIN, 1);
}

static int subframeIndex = 0;
void Thumby::WriteBufferGS(uint8_t* buffer, uint8_t* shaded, int bufferLength){
  absolute_time_t frame_time = get_absolute_time();
  // Write park commands, hold the display at row 1 until display ram has been written completly
  //sendCommand();
  //sendCommand();
  // With remaining time process some potential collisions
  // maximum cycles for n triangles is deterministic
  // process physics ticks at a lower framerate

  //frame_time = delayed_by_us(frame_time, FRAME_PERIOD);
  //sleep_until(frame_time);
}

// Send buffer to screen
void Thumby::writeBuffer(uint8_t* buffer, int bufferLength){
  sendCommand(SSD1306_COLUMNADDR);
  sendCommand(28);
  sendCommand(99);

  sendCommand(SSD1306_PAGEADDR);
  sendCommand(0x00);
  sendCommand(0x05);
  // TODO make this 4 drop the last page?

  digitalWrite(THUMBY_CS_PIN, 1);
  digitalWrite(THUMBY_DC_PIN, 1);
  digitalWrite(THUMBY_CS_PIN, 0);
  SPI.transfer(buffer, NULL, bufferLength);
  digitalWrite(THUMBY_CS_PIN, 1);
}


// Set the screen brightness to a value between 0 (off) and 127 (max)
void Thumby::setBrightness(uint8_t brightness){
  if(brightness>127){
    brightness = 127;
  }
  
  sendCommand(0x81);
  sendCommand(brightness);
}


// Return true if any buttons in the mask are detected as pressed
bool Thumby::isPressed(uint8_t mask){
  if(mask & (1 << 0) && !digitalRead(THUMBY_BTN_LDPAD_PIN)){
    return true;
  }else if(mask & (1 << 1) && !digitalRead(THUMBY_BTN_RDPAD_PIN)){
    return true;
  }else if(mask & (1 << 2) && !digitalRead(THUMBY_BTN_UDPAD_PIN)){
    return true;
  }else if(mask & (1 << 3) && !digitalRead(THUMBY_BTN_DDPAD_PIN)){
    return true;
  }else if(mask & (1 << 4) && !digitalRead(THUMBY_BTN_B_PIN)){
    return true;
  }else if(mask & (1 << 5) && !digitalRead(THUMBY_BTN_A_PIN)){
    return true;
  }

  return false;
}


// Pack dataBuf into packedBuf (adds 2 size bytes, 
// 1 checksum, and returns false if size too large to
// fit in packet, or too large to fit in packedBuf)
int8_t Thumby::linkPack(uint8_t* dataBuf, uint16_t dataBufLen, uint8_t* packedBuf, uint16_t packedBufLen){
  uint16_t packetLength = dataBufLen+3;

  // Check that the data length can be indexed by two bytes and
  // that it will fit into the packed buffer with header bytes
  if(dataBufLen > 512 || packetLength > packedBufLen){
    return -1;
  }

  // Prepare packet header
  packedBuf[0] = (dataBufLen >> 8) & 0xff;
  packedBuf[1] = dataBufLen & 0xff;
  packedBuf[2] = 0;

  // Generate checksum and copy data
  for(uint16_t b=0; b<dataBufLen; b++){
    packedBuf[2] ^= dataBuf[b];
    packedBuf[b+3] = dataBuf[b];
  }

  return packetLength;
}


// Unpack packedBuf to dataBuf (removes 1 checksum byte, and 2 size 
// bytes but returns false if checksum or size check fails, or if
// stripped packedBuf data cannot fit in dataBuf)
int8_t Thumby::linkUnpack(uint8_t* packedBuf, uint16_t packedBufLen, uint8_t* dataBuf, uint16_t dataBufLen){
  uint16_t dataLength = (packedBuf[0] << 8) + packedBuf[1];

  // Check that packet data will fit in data buffer and that
  // the received data length is the same as the actual received
  // packet length minus the 3 header bytes
  if(packedBufLen-3 > dataBufLen || dataLength != packedBufLen-3){
    return -1;
  }

  // Copy data and generate checksum off received data
  uint8_t checksum = 0;
  for(uint16_t b=0; b<dataLength; b++){
    dataBuf[b] = packedBuf[b+3];
    checksum ^= dataBuf[b];
  }

  // Return false if received and generated checksums are not the same
  if(checksum != packedBuf[2]){
    return -1;
  }

  return dataLength;
}


// Start playing a sound through the buzzer using pwm (does not block)
void Thumby::play(uint32_t freq, uint16_t duty){
  uint32_t wrap = clock_get_hz(clk_sys)/freq - 1;
  uint32_t level = (uint32_t)(wrap * (duty / 65535.0f));

  uint8_t slice_num = pwm_gpio_to_slice_num(THUMBY_AUDIO_PIN);
  pwm_set_wrap(slice_num, wrap);
  pwm_set_chan_level(slice_num, PWM_CHAN_A, level);
  
  // Set the PWM running
  pwm_set_enabled(slice_num, true);
}


// Turn off sound that's currently playing on pwm
void Thumby::stopPlay(){
  uint8_t slice_num = pwm_gpio_to_slice_num(THUMBY_AUDIO_PIN);
  
  // Set the PWM not running
  pwm_set_enabled(slice_num, false);
}