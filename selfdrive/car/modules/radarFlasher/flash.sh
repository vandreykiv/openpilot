#! /bin/sh
echo "================================="
echo "RADAR FLASH PROCESS STARTED"
echo "================================="
echo "NOTE: KEEP BRAKE PEDAL PRESSED UNTIL FLASH PROCESS IS COMPLETE"
echo "  "
echo "Starting the radar flash process..."
echo "  "
cd /data/openpilot/selfdrive/car/modules/radarFlasher
PYTHONPATH=/data/openpilot 
./patch.py --flash-firmware
ret=$?
if [ $ret -ne 0 ]; then
  echo "================================="
  echo "RADAR FLASH PROCESS FAILED"
  echo "================================="
  echo " Please check logs above for errors"
  echo " Please hit Reboot to return to OP"
  echo "An error occurred during flashing. Exiting..." >&2
  exit 1
fi
echo "================================="
echo "RADAR FLASH PROCESS COMPLETED"
echo "================================="
echo " Please hit Reboot to return to OP"
exit 0
