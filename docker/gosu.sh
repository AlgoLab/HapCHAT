#!/bin/bash

# Add local user
# with the same owner as /data
USER_ID=$(stat -c %u /data)
GROUP_ID=$(stat -c %g /data)

echo "Starting with UID:GID $USER_ID:$GROUP_ID"
groupadd -g "$GROUP_ID" group
useradd --shell /bin/bash -u "$USER_ID" -g group -o -c "" -m user
export HOME=/
chown --recursive "$USER_ID":"$GROUP_ID" /HapCHAT

exec gosu user "$@"
