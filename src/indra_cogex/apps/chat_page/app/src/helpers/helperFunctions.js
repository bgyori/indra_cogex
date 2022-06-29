/**
 * Helper functions for the chat page
 */

/* Get globally unique IDs */
let UUID = 0;
const getUUID = function () {
  return () => {
    UUID++;
    return UUID;
  };
};

export default {
  getUUID,
};
